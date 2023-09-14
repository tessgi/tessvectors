import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from matplotlib import cm

from lightkurve import TessQualityFlags
import numpy as np
import subprocess
import warnings
import logging
import requests
from bs4 import BeautifulSoup
import fitsio
import pandas as pd
import os

import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000000

import lightkurve as lk

import multiprocessing
from itertools import product

def get_ccd_centers(Sector, Camera, CCD):
    point=subprocess.run(["python","/Users/tapritc2/tessgi/tesspoint/tess-point/tess_stars2px.py",
                          "-r",str(Sector),str(Camera),str(CCD),
                          str(1024),str(1024)],capture_output=True,text=True)
    ra = float(point.stdout.split(' ')[0])
    dec = float(point.stdout.split(' ')[1])
    return ra, dec

#Get a timing array
def get_timing_midpoints_tpf(ra_c, dec_c, Sector, exptime):
    timing_benchmark=[]
    for ra, dec in zip(ra_c, dec_c):
        max_attempts=5
        while max_attempts:
            try: 
                if(exptime != 'ffi'):
                    results=lk.search_targetpixelfile(f"{ra} {dec}",radius=1e4,
                                                      sector=Sector,mission='TESS', 
                                                      exptime=exptime,limit=1, 
                                                      author='SPOC').download(quality_bitmask='none')
                else:
                    results=lk.search_tesscut(f"{ra} {dec}",
                                                      sector=Sector).download(cutout_size=(1,1),
                                                                              quality_bitmask='none')
                break
            except TimeoutError:
                    max_attempts -= 1
        
        #Get timecorr, correct for it
        #lightkurve PR?
        nonzero_mask = np.nonzero(results.time.value != 0)
        
        hdu=fits.open(results.path)
        timing_corr = hdu[1].data["TIMECORR"][nonzero_mask]
        timing_ccd=results.time.value[nonzero_mask] - timing_corr
        hdu.close
        
        
        if(len(timing_benchmark) == 0):
            timing_benchmark = timing_ccd
            cadences = results.cadenceno[nonzero_mask]#.value
            quality = results.quality[nonzero_mask]
        else:
            if(not np.allclose(timing_benchmark, timing_ccd, rtol=1e-6)):
                warnings.warn(f"Timing Cadence Mismatch in {exptime} TPF Data")
                Print(f"Length of Benchmark: {len(timing_benchmark)}")
                Print(f"Length of Comparison CCD: {len(timing_ccd)}")
                if(len(timing_benchmark) == len(timing_ccd)):
                    Print(f"Maximum Difference between the arrays: {max(abs(timing_benchmark-timing_ccd))}")
                #maybe set timing_benchmark = union of two timing sets here
                #swap above prints to logging-debug
    # Add in supplementrary Information
    TimingArr=np.zeros((4,len(timing_benchmark)), dtype=np.double)
    TimingArr[0] = timing_benchmark
    TimingArr[1] = timing_corr
    TimingArr[2] = cadences
    TimingArr[3] = quality
    return TimingArr
        
def get_camera_sector_cadences(Sector,Camera):
    cutoff_20s=27 # No 20s data before cycle_whatever
    CCD_List=[1]#(1,2,3,4)
    #write now we're checking that each CCD has the same timearr for each product
    #once we're done with that - we can can use only CCD 1 or whatever
    #IN PROGRESS
    
    ra_c=[]
    dec_c=[]
    #Get CCD Centers:
    for CCD in CCD_List:
        ra, dec = get_ccd_centers(Sector, Camera, CCD)
        ra_c.append(ra)
        dec_c.append(dec)
    
    midpoint_20s = None
    if(Sector >= cutoff_20s):
        #Get high-cadence 20s timing
        print(f"\t\tGetting 20s Cadence Midpoints")
        midpoint_20s = get_timing_midpoints_tpf(ra_c, dec_c, Sector, 'fast')
    #Get 120s TPF Timing
    midpoint_120s = get_timing_midpoints_tpf(ra_c, dec_c, Sector, 'short')
    print(f"\t\tGetting 120s Cadence Midpoints")
    #Get FFI Timing
    midpoint_ffi = get_timing_midpoints_tpf(ra_c, dec_c, Sector, 'ffi')
    print(f"\t\tGetting FFI Cadence Midpoints")

    return midpoint_20s, midpoint_120s, midpoint_ffi

def listFD(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def get_quat_data(Sector):
    base_url="https://archive.stsci.edu/missions/tess/engineering"
    filelist = [file for file in listFD(base_url,f"sector{Sector:02d}-quat.fits")]
    # Right now there is only 1 processing time per sector, this could change
    # For now assume this is a singular value
    # Later should potentially code some logic to use the laargest number/latest processing
    #QuatData=fitsio.FITS(filelist[0])
    print(f"\tGetting Quaternion Data for Sector: {Sector}")
    QuatData=fits.open(filelist[0])
    return QuatData

def get_emi_data(Sector):
    base_url="https://archive.stsci.edu/missions/tess/engineering"
    filelist = [file for file in listFD(base_url,f"sector{Sector:02d}-emi.fits")]
    # Right now there is only 1 processing time per sector, this could change
    print(f"\tGetting EMI Data for Sector: {Sector}")
    EMIData=fits.open(filelist[0])
    return EMIData

def crm_bin(quats):
    #we could re-write this to be O(N) if we just step through the array checking highest/lowest values
    if (len(quats) > 2):
        return np.median(np.sort(quats)[1:len(quats)-1])
    else:
        return

def Bin_Quat_Cadence(cadence, CameraData, Camera):
    ExpTime = np.double(abs(cadence[1]-cadence[0])) # Should this be hard-coded instead? prbly med
    #print(ExpTime)
    SubArr_Len = int(2 * (ExpTime * (60. * 60. * 24.) / 2.))
    Time_Len = len(cadence)
    Source_Len = len(CameraData['Time'])
    BinnedQuats=np.zeros((21, Time_Len), dtype=np.double) #( 8 Quat Vectors * 3 + 2 (Tmin/max) + #GS + MSTot + FoM)
    submask_length = 2 * ExpTime
    min_ind = int(0)
    max_ind = SubArr_Len # Days -> Seconds / 2 S per Quat Exp
    Flag_ZeroMask = 0
    #print(f"Index Time Length: {Time_Len}")
    #print(f"Source Time Length: {len(CameraData)}")

    
    for i in range(Time_Len - 1) :
        #print(f"min: {min_ind} max: {max_ind} max_len: {Time_Len} SubArr_Len: {SubArr_Len}")
        #print(f"i: {i} ")
        
        if(i < (Time_Len - 1)):
            while CameraData['Time'][max_ind] < cadence[i+1]: 
                max_ind = max_ind + SubArr_Len
        if(i > 1):
            while CameraData['Time'][min_ind]  > cadence[i-1]: 
                min_ind = min_ind - 1        

        SubArr = CameraData[min_ind:max_ind]

        mask = np.nonzero(abs(np.double(SubArr['Time'] - cadence[i])) < np.double(np.double(0.5) * ExpTime))
        
        if(len(mask[0]) == 0):
            Flag_ZeroMask = Flag_ZeroMask + 1
            #print(f"\t\t\t\tWarning: 0 Quaternions found in the time bin!")
            #print(f"Min: {CameraData['Time'][min_ind]}")
            #print(f"Max: {CameraData['Time'][max_ind]}")
            #print(f"Cadence: {cadence[i]}")
            #print(f"ExpTime: {ExpTime}")
            #print(f"ArrEq: {abs(np.double(SubArr['Time'] - cadence[i])) < np.double(np.double(0.5) * ExpTime)}")
            #print(f"DMin: {cadence[i] - 0.5 * ExpTime}")
            #print(f"DMax: {cadence[i] + 0.5 * ExpTime}")
            #print(f"SubArr: {SubArr['Time']}")
            #print(f"Abs(SubArr-Cadence): {abs(np.double(SubArr['Time'] - cadence[i]))}")

        #print(SubArr['Time'] - cadence[i])
        #print(cadence[i])
        #print(mask)
        #print(f"Mask: {mask}")
        #print(f"Mask_Len: {len(mask[0])}")
        if(len(mask[0]) != 0):
            BinnedQuats[0][i] = np.min(SubArr['Time'][mask]) #Min Midpoint Binned
            BinnedQuats[1][i] = np.max(SubArr['Time'][mask]) #Max Midpoint Binned
            BinnedQuats[2][i] = np.min(SubArr[f"C{Camera:d}_FOM"][mask]) # Figure of Merit 
            BinnedQuats[3][i] = np.min(SubArr[f"C{Camera:d}_NUM_GSUSED"][mask]) # Number of Guide Stars Used # remove this or FoM
            BinnedQuats[4][i] = len(mask[0]) # N quats used - diagnose wierdness post binning
            
            BinnedQuats[5][i] =   np.median(SubArr[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[6][i] =   np.std(SubArr[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[7][i] =   sigma_clipped_stats(SubArr[f"C{Camera:d}_Q1"][mask])[2]
            BinnedQuats[8][i] =   crm_bin(SubArr[f"C{Camera:d}_Q1"][mask])

            BinnedQuats[9][i] =   np.median(SubArr[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[10][i] =  np.std(SubArr[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[11][i] =  sigma_clipped_stats(SubArr[f"C{Camera:d}_Q2"][mask])[2]
            BinnedQuats[12][i] =  crm_bin(SubArr[f"C{Camera:d}_Q2"][mask])
            
            BinnedQuats[13][i] =  np.median(SubArr[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[14][i] =  np.std(SubArr[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[15][i] =  sigma_clipped_stats(SubArr[f"C{Camera:d}_Q3"][mask])[2]
            BinnedQuats[16][i] =  crm_bin(SubArr[f"C{Camera:d}_Q3"][mask])
            
            BinnedQuats[17][i] =  np.median(SubArr[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[18][i] =  np.std(SubArr[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[19][i] =  sigma_clipped_stats(SubArr[f"C{Camera:d}_Q4"][mask])[2]
            BinnedQuats[20][i] =  crm_bin(SubArr[f"C{Camera:d}_Q4"][mask])
           
            min_ind = min_ind + int(np.max(mask[0]) + 1)

        else:
            for j in range(0,16):
                BinnedQuats[j][i]=np.NaN
                #Sector 43 - Camera 2 not being used, only 4
                #Running into issues where there are no quaternions
                #Substitute Nan's?
                #Include data quality flags, check to see if this is in known bad data and
                #throw a harder error?
        # Add Cadence #
        # Quality Masks
        #print(max(mask), SubArr_Len)
        max_ind = int(min_ind)+int(SubArr_Len)
        if(max_ind > Source_Len): 
            max_ind=Source_Len-1
        
    if(Flag_ZeroMask > 0): 
        print(f"\t\t\t\tWARNING: Zero Quaternions found for {Flag_ZeroMask} observed cadences")
    return BinnedQuats

def EMI_Extension_Dict(ext_string):
    # This could be done algorithmically, but thought it was less error prone this way
    # I broke this out of the function b/c it was big & ugly
     emi_dict={
       "emi.earth_distance": 1,
       "emi.earth_sc_elevation": 2,
       "emi.earth_sc_azimuth": 3,
       "emi.moon_distance": 4,
       "emi.moon_sc_elevation": 5,
       "emi.moon_sc_azimuth": 6,
       "emi.cam1_earth_boresight_angle": 7,
       "emi.cam1_earth_azimuth": 8,
       "emi.cam1_moon_boresight_angle": 9,
       "emi.cam1_moon_azimuth": 10,
       "emi.cam2_earth_boresight_angle": 11,
       "emi.cam2_earth_azimuth": 12,
       "emi.cam2_moon_boresight_angle": 13,
       "emi.cam2_moon_azimuth": 14,
       "emi.cam3_earth_boresight_angle": 15,
       "emi.cam3_earth_azimuth": 16,
       "emi.cam3_moon_boresight_angle": 17,
       "emi.cam3_moon_azimuth": 18,
       "emi.cam4_earth_boresight_angle": 19,
       "emi.cam4_earth_azimuth": 20,
       "emi.cam4_moon_boresight_angle": 21,
       "emi.cam4_moon_azimuth": 22
       }
     return emi_dict[ext_string]

def Bin_EMI_Cadence(cadence, EMIData, Camera,type):
    typedict = {1:'020', 2:'120', 3:'FFI'}
    ExpTime = np.double(abs(cadence[1]-cadence[0])) # Should this be hard-coded instead? prbly med
    if (isinstance(type, int)):
        type=typedict[type]
    match type:
        case "020":
            SubArr_Len=int(4)
        case "120":
            SubArr_Len=int(4)
        case "FFI":
            SubArr_Len= int((ExpTime / (EMIData[1].data['TIME'][1] - EMIData[1].data['TIME'][0])) * 4)
        case _:
            print("Error - Bin_EMI_Cadence can't match cadence type")
            #throw a logging warning instead
    
    Time_Len = len(cadence)
    Source_Len = len(EMIData[1].data['Time'])
    BinnedEMI=np.zeros((6, Time_Len), dtype=np.double) #( (Earth + Moon) ** (Distance + Elev. + Azim.)
    min_ind = int(0)
    max_ind = int(SubArr_Len)
   
    for i in range(Time_Len - 1) :
        if(i < (Time_Len - 2)):
            while EMIData[1].data['Time'][max_ind] < cadence[i+1]: 
                max_ind = max_ind + SubArr_Len

        if(i >= 1):
            while EMIData[1].data['Time'][min_ind]  > cadence[i-1]: 
                min_ind = min_ind - 1 
       

        SubArr = EMIData[1].data['Time'][min_ind:max_ind]
        match_ind = np.argmin(abs(np.double(SubArr - cadence[i])))
        #print(f"i: {i} min_ind: {min_ind} max_ind:{max_ind} match_ind:{match_ind} Time_Len:{Time_Len} Subarr_Len:{SubArr_Len} ")

        if(isinstance(min_ind, int)):
            # Earth Distance
            BinnedEMI[0][i] = EMIData[EMI_Extension_Dict('emi.earth_distance')].data[min_ind:max_ind]['VALUE'][match_ind]
            # Earth Camera Elev
            BinnedEMI[1][i] = EMIData[EMI_Extension_Dict(f'emi.cam{Camera}_earth_boresight_angle')].data[min_ind:max_ind]['VALUE'][match_ind] 
            # EarthCamera Azim
            BinnedEMI[2][i] = EMIData[EMI_Extension_Dict(f'emi.cam{Camera}_earth_azimuth')].data[min_ind:max_ind]['VALUE'][match_ind]
            # Moon Distance
            BinnedEMI[3][i] = EMIData[EMI_Extension_Dict('emi.moon_distance')].data[min_ind:max_ind]['VALUE'][match_ind]
            # Moon Camera Elev
            BinnedEMI[4][i] = EMIData[EMI_Extension_Dict(f'emi.cam{Camera}_moon_boresight_angle')].data[min_ind:max_ind]['VALUE'][match_ind] 
            # Moon Camera Azim
            BinnedEMI[5][i] = EMIData[EMI_Extension_Dict(f'emi.cam{Camera}_moon_azimuth')].data[min_ind:max_ind]['VALUE'][match_ind]
        else:
            print("Bin_EMI_Cadence can't match a datapoint")
            #throw a real logging.warning
        
        min_ind = min_ind + int(match_ind)
        max_ind = int(min_ind) + int(SubArr_Len)
        if(max_ind > Source_Len): 
            max_ind=Source_Len-1
 
    return BinnedEMI

def Bin_Quat_Camera(CameraData,TimeArr, Camera):

    Bin20 = None   
    if(TimeArr[0] is not None):
        print("\t\tBinning 20s Quaternion Data")
        Bin20 = Bin_Quat_Cadence(TimeArr[0][0], CameraData, Camera)
    print("\t\tBinning 120s Quaternion Data")
    Bin120 = Bin_Quat_Cadence(TimeArr[1][0], CameraData, Camera)
    print("\t\tBinning FFI Quaternion Data")
    BinFFI = Bin_Quat_Cadence(TimeArr[2][0], CameraData, Camera)
        
    return Bin20, Bin120, BinFFI

def Bin_EMI_Camera(EMIData,TimeArr, Camera):
    EMI20 = None   
    if(TimeArr[0] is not None):
        print("\t\tBinning 20s EMI Data")
        EMI20 = Bin_EMI_Cadence(TimeArr[0][0], EMIData, Camera,'020')
    print("\t\tBinning 120s EMI Data")
    EMI120 = Bin_EMI_Cadence(TimeArr[1][0], EMIData, Camera, '120')
    print("\t\tBinning FFI EMI Data")
    EMIFFI = Bin_EMI_Cadence(TimeArr[2][0], EMIData, Camera, 'FFI')
        
    return EMI20, EMI120, EMIFFI

def write_vector_sector_camera(BinnedQuats, BinnedEMI, TimeArr, Sector, Camera):
    typedict = {1:'020', 2:'120', 3:'FFI'}
    Binned_Dir='TESSVectors_products'

    for Quat, EMI, Time, i in zip(BinnedQuats, BinnedEMI, TimeArr, [1,2,3]):
        from datetime import date

        fname = f"{Binned_Dir}/{typedict[i]}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{typedict[i]}.csv"
        print(f"\t\tWriting Sector: {Sector} Camera: {Camera} Cadence {typedict[i]} to:")
        print(f"\t\t\t to {fname}")
        
        if (Time is not None):
            # First write the File Header
            # this will over-write by default
            with open(fname, 'w') as f:
                f.write(f"# TESS Quaternions downsampled to end-user cadences\n")
                f.write(f"# Sector: \t{Sector}\n")
                f.write(f"# Camera: \t{Camera}\n")
                f.write(f"# Cadence:\t{typedict[i]}\n")
                f.write(f"# This file contains TESS quaternions that have been downsampled from the TESS\n")
                f.write(f"# native 2-second exposures to the end-user product cadence (e.g .20s/120s/FFI)\n")
                f.write(f"# For more information See:\n")
                f.write(f"#     - The github repo that created this file at https://github.com/tylerapritchard/TESSQuats\n")# This will presumably move to tessgi at some point
                f.write(f"#     - The TESS Instrument Handbook at https://archive.stsci.edu/missions-and-data/tess\n")
                f.write(f"#     - The original quaternion engineering files at https://archive.stsci.edu/missions-and-data/tess/data-products\n")
                f.write(f"# Please check the TESS Sector Data Release Notes at https://archive.stsci.edu/tess/tess_drn.html\n")
                f.write(f"# for information regarding the telescope pointing and tracking (in addition to other image quality information) as\n")
                f.write(f"# for some sectors/cameras alternative camera quaternions are substituted for a variety of programattic reasons.\n\n")
                f.write(f"# Column Descriptions:\n")
                f.write(f"# Cadence #: Cadence index from the source tpf\n")
                # NOTE TO SELF WE SHOULD PROBABLY RE-CONVERT THE TIMES TO BE LIKE LK TIMES FOR THE END USER
                f.write(f"# MidTime: The exposure midpoint in spacecraft time (i.e. tpf.time - tpf.timecorr)\n")
                f.write(f"# TimeCorr: The Time Correction for spacecraft to Barycentric time at that cadence\n")
                f.write(f"# Quality: Quality bitmask at that cadence from the tpf - see the TESS Instrument handbook for more info\n")
                # We could remove ExpTime/Sector/Camera and just keep them in the metadata, but I could see them being
                # useful like this
                f.write(f"# ExpTime: The final cadence binning (20s/120s/FFI)\n")
                f.write(f"# Sector: The TESS observing Sector for the source data\n")
                f.write(f"# Camera: The TESS camera for the source data\n")
                f.write(f"# Quat_Start: The timestamp of the earliest quaternion used in the bin\n")
                f.write(f"# Quat_Stop: The timestamp of the last quaternion used in the bin\n")
                f.write(f"# Quat_MIN_FOM: The worst Figure of Merit from the source quaternions\n")
                f.write(f"# Quat_MIN_NUM_GSUSED: The lowest number of guide stars used in the source quaternions\n")
                f.write(f"# Quat_NBinned: The number of quaternions binned into this final result.\n")
                f.write(f"# Quat[1-4]_Med: The Quaternion #[1-4] median value from the binned values \n")
                f.write(f"# Quat[1-4]_StdDev: The standard deviation of Quaternion #[1-4] binned values\n")
                f.write(f"# Quat[1-4]_SigClip: The Sigma-Clipped Standard Deviation of Quaternion #[1-4] binned values\n")
                f.write(f"# Quat[1-4]_CRM_Med: The Quaternion #[1-4] median value with the highest and lowest values excluded \n\n")
                f.write(f"# Earth_Distance: Distance to Earth in Re \n\n")            
                f.write(f"# Earth_Camera_Angle: Angle of Earth from Camera Boresight in Degrees \n\n")            
                f.write(f"# Earth_Camera_Azimuth: Azimuth of Earth around Camera Boresight in Degrees \n\n")            
                f.write(f"# Moon_Distance: Distance to Moon in Re \n\n")            
                f.write(f"# Moon_Camera_Angle: Angle of Moon from Camera Boresight in Degrees \n\n")            
                f.write(f"# Moon_Camera_Azimuth: Azimuth of Moon around Camera Boresight in Degrees \n\n")            

                f.write(f"# Processing Date-Time: {date.today()}\n\n")

            df = pd.DataFrame(data={'Cadence':Time[2], 'MidTime':Time[0] + Time[1], 'TimeCorr':Time[1],
                                    'Quality':Time[3], 'ExpTime':[typedict[i]] * len(Time[0]), 
                                    'Sector': [Sector] * len(Time[0]), 'Camera': Camera  * len(Time[0]),
                                    'Quat_Start':Quat[0], 'Quat_Stop':Quat[1], 
                                    'Quat_MIN_FOM':Quat[2],'Quat_MIN_NUM_GSUSED':Quat[3],'Quat_NBinned': Quat[4], 
                                    'Quat1_Med': Quat[5], 'Quat1_StdDev': Quat[6], 'Quat1_StdDev_SigClip':Quat[7],'Quat1_CRM_Med': Quat[8], 
                                    'Quat2_Med': Quat[9], 'Quat2_StdDev': Quat[10],'Quat2_StdDev_SigClip':Quat[11],'Quat2_CRM_Med': Quat[12], 
                                    'Quat3_Med': Quat[13],'Quat3_StdDev': Quat[14],'Quat3_StdDev_SigClip':Quat[15],'Quat3_CRM_Med': Quat[16], 
                                    'Quat4_Med': Quat[17],'Quat4_StdDev': Quat[18],'Quat4_StdDev_SigClip':Quat[19],'Quat4_CRM_Med': Quat[20],
                                    'Earth_Distance': EMI[0], 'Earth_Camera_Angle': EMI[1], 'Earth_Camera_Azimuth': EMI[2],
                                    'Moon_Distance': EMI[3], 'Moon_Camera_Angle': EMI[4], 'Moon_Camera_Azimuth': EMI[5] 
                                    }
                             )
            df.astype({'Cadence':int, 
                       'Quality':int, 
                       'Sector': int, 
                       'Camera':int}
                       #'Quat_MIN_NUM_GSUSED':int,
                       #'Quat_NBinned': int}
                    ).to_csv(fname, index=False, mode = "a")
    
def bin_Sector(Sector):
    print(f"Starting Sector: {Sector}")
    QuatData = get_quat_data(Sector)
    EMIData = get_emi_data(Sector)
    for Camera in [1,2,3,4]:
        print(f"\tStarting Camera: {Camera}")
        TimeArr = get_camera_sector_cadences(Sector, Camera)
        BinnedQuats = Bin_Quat_Camera(QuatData[Camera].data, TimeArr, Camera)
        BinnedEMI = Bin_EMI_Camera(EMIData, TimeArr, Camera)
        write_vector_sector_camera(BinnedQuats, BinnedEMI, TimeArr, Sector, Camera)

def quality_to_color(qual):
    
    default = TessQualityFlags.create_quality_mask(qual, bitmask='default')
    hard = TessQualityFlags.create_quality_mask(qual, bitmask='hard')
    #hardest = TessQualityFlags.create_quality_mask(qual, bitmask='hardest')

    carr = np.zeros(len(qual))
    #carr[~ hardest] = 0.33
    carr[~ hard] = 0.5
    carr[~ default] = 1
    
    return carr


def plot_quat(axs, time, quat, dev, qual,QuatLabel):
        
        nsigma = 1
        
        #norm = np.nanmedian(np.double(quat[~np.isfinite(quat)]))
        norm=np.nanmedian(quat)
        norm_quat = quat / norm
        norm_dev = dev / norm
        #mean_norm_dev = np.nanmedian(dev[~np.isfinite(dev)] / norm)
        mean_norm_dev = np.nanmedian(dev / norm)
        tmin=np.min(time)
        tmax=np.max(time)
        
        ymin = np.nanmedian(norm_quat) - 3 * mean_norm_dev
        ymax = np.nanmedian(norm_quat) + 3 * mean_norm_dev
        
        im = axs.imshow(np.vstack((qual,)), 
                        extent=(tmin, tmax, ymin, ymax), 
                        interpolation='nearest', aspect='auto', 
                        cmap=cm.PuRd, vmax=1)  
        
        axs.scatter(time, norm_quat, s=3, 
                    label = 'Median Quaternion Value', color='k')
        
        axs.fill_between(time, 
                            norm_quat - (nsigma * norm_dev),
                            norm_quat + (nsigma * norm_dev),
                            alpha = 0.6, color='grey')
        
        axs.set_xlim(tmin, tmax)        
        axs.set_ylim(ymin, ymax)
        
        axs.tick_params(axis='x', labelsize=0)
        
        axs.set_ylabel(f"{QuatLabel}", weight='bold', size=24 )
        axs.set_yticks([])
        
        return im

def create_diagnostic_timeseries(Sector, Camera, Cadence):

    typedict = {1:'020', 2:'120', 3:'FFI'}
    if(type(Cadence) != str):
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    Binned_Dir='TESSVectors_products'
    fname = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"
    
    nplots=3
    if(os.path.isfile(fname)):
        quatdf = pd.read_csv(fname, comment='#',index_col=False)
        qual_im = quality_to_color(quatdf.Quality)

        fig, axs = plt.subplots(nplots,1,figsize=(15,nplots*10))
        im = plot_quat(axs[0], quatdf.MidTime, quatdf.Quat1_Med, quatdf.Quat1_StdDev, qual_im, 'Quaternion 1')
        im = plot_quat(axs[1], quatdf.MidTime, quatdf.Quat2_Med, quatdf.Quat2_StdDev, qual_im, 'Quaternion 2')
        im = plot_quat(axs[2], quatdf.MidTime, quatdf.Quat3_Med, quatdf.Quat3_StdDev, qual_im, 'Quaternion 3')
        plt.subplots_adjust(hspace=0)
        axs[0].set_title(f"TESS Sector {Sector} Camera {Camera} Quaternions",
                         weight="bold", size=26)
        axs[-1].tick_params(axis='x', labelsize=18)
        axs[-1].set_xlabel("TESS BTJD", weight='bold', size=24)
        
        cax = plt.axes([0.92, 0.11, 0.075, 0.77])
        cbar = plt.colorbar(mappable = im, cax=cax, 
                            ticks=[0,0.5,1])
        cbar.ax.set_yticklabels(['Unflagged', 'Aggressive', 'Conservative'], 
                                size=18)
        cbar.set_label("Data Flagging Level (Lower Is Better)", 
                       size=24, weight='bold')
        
        #plt.tight_layout()
        fout = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_Quat.png"
        plt.show()
        plt.savefig(fout, dpi=300,bbox_inches='tight', rasterize=True)
        plt.close(fig)

def plot_lsperiodogram(ax, time, median, std, QuatLabel):
    lc = lk.LightCurve(data={'time': time , 'flux': std})
    ls=lc.to_periodogram(maximum_period=13.5)
    ls.plot(ax=ax,lw=0.1, color='k', ylabel=" ") 

    ax.set_ylabel(f"{QuatLabel} Power", weight='bold', size=24 )
    #ax.set_yticks([])
    ax.tick_params(axis='y', labelsize=18)
    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_ylim(min(ls.power),max(ls.power))
    
    return

def create_diagnostic_periodogram(Sector, Camera, Cadence):
    
    typedict = typedict = {1:'020', 2:'120', 3:'FFI'}
    if(type(Cadence) != str):
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    Binned_Dir='TESSVectors_products'
    fname = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"
    
    nplots=3
    if(os.path.isfile(fname)):
        quatdf = pd.read_csv(fname, comment='#',index_col=False)   
        fig, axs = plt.subplots(nplots,1,figsize=(15,nplots*10))
        plot_lsperiodogram(axs[0], quatdf.MidTime, quatdf.Quat1_Med, quatdf.Quat1_StdDev, 'Quaternion 1')
        plot_lsperiodogram(axs[1], quatdf.MidTime, quatdf.Quat2_Med, quatdf.Quat2_StdDev, 'Quaternion 2')
        plot_lsperiodogram(axs[2], quatdf.MidTime, quatdf.Quat3_Med, quatdf.Quat3_StdDev, 'Quaternion 3')
        plt.subplots_adjust(hspace=0)

        axs[0].set_title(f"TESS Sector {Sector} Camera {Camera} Quaternion Power Spectra",
                         weight="bold", size=26)
        axs[-1].tick_params(axis='x', labelsize=18)
        axs[-1].set_xlabel("Period [days]", weight='bold', size=24)

        forward = lambda x: x * 24. * 3600.
        inverse = lambda x: x / 24. / 3600.
        ax2 = axs[0].secondary_xaxis('top', functions=(forward, inverse))
        ax2.set_xlabel("Period [seconds]", weight='bold', size=24)
        ax2.tick_params(axis='x', labelsize=18)

        fout = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_QuatPower.png"
        plt.savefig(fout, dpi=300,bbox_inches='tight')
        plt.close(fig)

def create_diagnostic_emi(Sector, Camera, Cadence):
    import matplotlib as mpl
    mpl.rcParams['agg.path.chunksize'] = 10000000


    typedict = typedict = {1:'020', 2:'120', 3:'FFI'}
    if(type(Cadence) != str):
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    Binned_Dir='TESSVectors_products'
    fname = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"
    
    nplots=3
    if(os.path.isfile(fname)):
        quatdf = pd.read_csv(fname, comment='#',index_col=False)

        qual_im = quality_to_color(quatdf.Quality)   
        fig, axs = plt.subplots(nplots,1,figsize=(15,nplots*10))
        
        axs[0].scatter(quatdf.MidTime,quatdf.Earth_Distance, 
                       label = "Earth", color='seagreen')
        axs[0].scatter(quatdf.MidTime,quatdf.Moon_Distance,  
                       label = "Moon",  color='darkturquoise')
        axs[0].legend(prop={'weight':'bold', 'size': 24},
                      scatterpoints=3, markerscale=2)

        tmin=np.min(quatdf.MidTime)
        tmax=np.max(quatdf.MidTime)

        ymin=min(min(quatdf.Earth_Distance), min(quatdf.Moon_Distance))
        ymax=max(max(quatdf.Earth_Distance), max(quatdf.Moon_Distance))
        im = axs[0].imshow(np.vstack((qual_im,)), 
                           extent=(tmin, tmax, ymin, ymax), 
                           interpolation='nearest', aspect='auto', 
                           cmap=cm.PuRd, vmax=1)  

        axs[1].scatter(quatdf.MidTime,quatdf.Earth_Camera_Angle, 
                       label = "Earth", color='seagreen')
        axs[1].scatter(quatdf.MidTime,quatdf.Moon_Camera_Angle, 
                       label = "Moon",   color='darkturquoise')   
        axs[1].legend(prop={'weight':'bold', 'size': 24},
                      scatterpoints=3, markerscale=2)

        ymin=min(min(quatdf.Earth_Camera_Angle), min(quatdf.Moon_Camera_Angle))
        ymax=max(max(quatdf.Earth_Camera_Angle), max(quatdf.Moon_Camera_Angle))
        im = axs[1].imshow(np.vstack((qual_im,)), 
                           extent=(tmin, tmax, ymin, ymax), 
                           interpolation='nearest', aspect='auto', 
                           cmap=cm.PuRd, vmax=1)  


        axs[2].scatter(quatdf.MidTime,quatdf.Earth_Camera_Azimuth, 
                       label = "Earth", color='seagreen')
        axs[2].scatter(quatdf.MidTime,quatdf.Moon_Camera_Azimuth, 
                       label = "Moon",   color='darkturquoise')   
        axs[2].legend(prop={'weight':'bold', 'size': 24},
                      scatterpoints=3, markerscale=2)

        ymin=min(min(quatdf.Earth_Camera_Azimuth), min(quatdf.Moon_Camera_Azimuth))
        ymax=max(max(quatdf.Earth_Camera_Azimuth), max(quatdf.Moon_Camera_Azimuth))
        im = axs[2].imshow(np.vstack((qual_im,)), 
                           extent=(tmin, tmax, ymin, ymax), 
                           interpolation='nearest', aspect='auto', 
                           cmap=cm.PuRd, vmax=1)  

        axs[0].set_ylabel("Distance", weight='bold', size=24)
        axs[1].set_ylabel("Camera {Camera} Angle", weight='bold', size=24)
        axs[2].set_ylabel("Camera {Camera} Azimuth", weight='bold', size=24)

        axs[2].set_xlabel("TESS BTJD", weight='bold', size=24)
        
        axs[0].tick_params(axis='y', labelsize=18)
        axs[1].tick_params(axis='y', labelsize=18)
        axs[2].tick_params(axis='y', labelsize=18)

        axs[2].tick_params(axis='x', labelsize=18)
        plt.subplots_adjust(hspace=0)

        cax = plt.axes([0.92, 0.11, 0.075, 0.77])
        cbar = plt.colorbar(mappable = im, cax=cax, 
                            ticks=[0,0.5,1])
        cbar.ax.set_yticklabels(['Unflagged', 'Aggressive', 'Conservative'], 
                                size=18)
        cbar.set_label("Data Flagging Level (Lower Is Better)", 
                       size=24, weight='bold')
        
        fout = f"{Binned_Dir}/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_emi.png"
        plt.savefig(fout, dpi=300,bbox_inches='tight')
        #plt.show()
        plt.close(fig)

def Create_Diagnostics_Sector(Sector):
    #CameraCadence
    #Sector, Camera, Cadence = SectorCameraCadence
    print(f"Creating Diagnostics for Sector: {Sector}")
    for Camera in [1,2,3,4]:
        for Cadence in [1,2,3]:
            create_diagnostic_timeseries(Sector, Camera, Cadence)
            #Should I create the periodograms from the "raw" 2s data?  probably?
            create_diagnostic_periodogram(Sector, Camera, Cadence)
            create_diagnostic_ephemerides(Sector, Camera, Cadence)

def TESSVectors_process_sector(Sector):
    bin_Sector(Sector)
    try:
        Create_Diagnostics_Sector(Sector)
    except:
        print("\t\t\t Warning, Plotting failed") # add a real warning
        pass

    return(Sector, True)

def run_bulk_diagnostics():
    sector_list = range(1,65)
    camera_list = range(1,4)
    cadence_list = range(1,3)
    inlist=list(product(sector_list,camera_list,cadence_list))
    
    from multiprocessing.pool import Pool
    from itertools import product

    pool=Pool(processes=7)
    res=[]
    for result in pool.map(create_diagnostics_bulk, inlist):
        res=[res,result]
    pool.close()

def run_bulk_quats(processes=7):
    if (not processes):
        processes = 7
    sector_list = range(1,65)

    from multiprocessing.pool import Pool

    pool=Pool(processes=processes)
    res=[]
    for result in pool.map(bin_Sector, sector_list):
        res=[res,result]
    pool.close()

def run_bulk_processing():
    from multiprocessing.pool import Pool
    import multiprocessing as mp
    sector_list = range(1,65)

    pool=Pool(processes=mp.cpu_count())
    res=[]
    for result in pool.map(TESSVectors_process_sector, sector_list):
        res=[res,result]
    pool.close()