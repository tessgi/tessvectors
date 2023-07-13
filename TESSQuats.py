#Camera, Sector = 1, 1
#Processing steps: grab a 20s, 120s, FFI cutout to get midpoints

#Get the center of a CCD for a given Camera/Sector/CCD
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
import subprocess
import warnings
import logging
import astropy.io.fits as fits
import requests
from bs4 import BeautifulSoup
import fitsio
import pandas as pd
from astropy.stats import sigma_clipped_stats

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

def Bin_Cadence(cadence, CameraData, Camera):
    ExpTime = np.double(abs(cadence[1]-cadence[0])) # Should this be hard-coded instead? prbly med
    #print(ExpTime)
    SubArr_Len = int(2 * (ExpTime * (60. * 60. * 24.) / 2.))
    Time_Len = len(cadence)
    Source_Len = len(CameraData['Time'])
    BinnedQuats=np.zeros((17, Time_Len), dtype=np.double) #( 8 Quat Vectors * 3 + 2 (Tmin/max) + #GS + MSTot + FoM)
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
            BinnedQuats[0][i] = np.min(CameraData['Time'][mask]) #Min Midpoint Binned
            BinnedQuats[1][i] = np.max(CameraData['Time'][mask]) #Max Midpoint Binned
            BinnedQuats[2][i] = np.min(CameraData[f"C{Camera:d}_FOM"][mask]) # Figure of Merit 
            BinnedQuats[3][i] = np.min(CameraData[f"C{Camera:d}_NUM_GSUSED"][mask]) # Number of Guide Stars Used # remove this or FoM
            BinnedQuats[4][i] = len(mask[0]) # N quats used - diagnose wierdness post binning
            
            BinnedQuats[5][i] =   np.median(CameraData[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[6][i] =   np.std(CameraData[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[7][i] =   sigma_clipped_stats(CameraData[f"C{Camera:d}_Q1"][mask])[2]
            
            BinnedQuats[8][i] =   np.median(CameraData[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[9][i] =   np.std(CameraData[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[10][i] =  sigma_clipped_stats(CameraData[f"C{Camera:d}_Q2"][mask])[2]
            
            BinnedQuats[11][i] =  np.median(CameraData[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[12][i] =  np.std(CameraData[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[13][i] =  sigma_clipped_stats(CameraData[f"C{Camera:d}_Q3"][mask])[2]
            
            BinnedQuats[14][i] =  np.median(CameraData[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[15][i] =  np.std(CameraData[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[16][i] =  sigma_clipped_stats(CameraData[f"C{Camera:d}_Q4"][mask])[2]
           
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
        max_ind = int(min_ind+SubArr_Len)
        if(max_ind > Source_Len): 
            max_ind=Source_Len-1
        
    if(Flag_ZeroMask > 0): 
        print(f"\t\t\t\tWARNING: Zero Quaternions found for {Flag_ZeroMask} observed cadences")
    return BinnedQuats

def Bin_Camera(CameraData,TimeArr, Camera):

    Bin20 = None   
    if(TimeArr[0] is not None):
        print("\t\tBinning 20s Data")
        Bin20 = Bin_Cadence(TimeArr[0][0], CameraData, Camera)
    print("\t\tBinning 120s Data")
    Bin120 = Bin_Cadence(TimeArr[1][0], CameraData, Camera)
    print("\t\tBinning FFI Data")
    BinFFI = Bin_Cadence(TimeArr[2][0], CameraData, Camera)
        
    return Bin20, Bin120, BinFFI

def write_quat_sector_camera(BinnedQuats, TimeArr, Sector, Camera):
    Binned_Dir='binned_quats'
    typedict = {1:'fast', 2:'short', 3:'long'}

    for Quat, Time, i in zip(BinnedQuats, TimeArr, [1,2,3]):
        fname = f"{Binned_Dir}/tess_BinQuats_S{Sector:03d}_C{Camera}_{typedict[i]}.csv"
        print(f"\t\tWriting Sector: {Sector} Camera: {Camera} Cadence {typedict[i]} to:")
        print(f"\t\t\t to {fname}")
        
        if (Time is not None):
            df = pd.DataFrame(data={'MidTime':Time[0], 'TimeCorr':Time[1],
                                    'CadenceNo':Time[2], 'Quality':Time[3], 'ExpTime':[typedict[i]] * len(Time[0]), 
                                    'Sector': [Sector] * len(Time[0]), 'Camera': Camera  * len(Time[0]),
                                   'Quat_Start':Quat[0], 'Quat_Stop':Quat[1], 
                                   'Quat_MIN_FOM':Quat[2],'Quat_MIN_NUM_GSUSED':Quat[3],
                                   'Quat_NBinned': Quat[4], 
                                   'Quat1_Med': Quat[5], 'Quat1_StdDev': Quat[6], 'Quat1_SigClip':Quat[7],
                                   'Quat2_Med': Quat[8], 'Quat2_StdDev': Quat[9], 'Quat2_SigClip':Quat[10], 
                                   'Quat3_Med': Quat[11],'Quat3_StdDev': Quat[12],'Quat3_SigClip':Quat[13],
                                   'Quat4_Med': Quat[14],'Quat4_StdDev': Quat[15],'Quat4_SigClip':Quat[16]}
                             )
            df.to_csv(fname)
    
def bin_Sector(Sector):
    print(f"Starting Sector: {Sector}")
    QuatData = get_quat_data(Sector)
    for Camera in [1,2,3,4]:
        print(f"\tStarting Camera: {Camera}")
        TimeArr = get_camera_sector_cadences(Sector, Camera)
        BinnedQuats = Bin_Camera(QuatData[Camera].data, TimeArr, Camera)
        write_quat_sector_camera(BinnedQuats, TimeArr, Sector, Camera)
        return Sector
