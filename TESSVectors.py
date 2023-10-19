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
import pandas as pd
import os
import sys
import matplotlib as mpl

mpl.rcParams["agg.path.chunksize"] = 10000000

import lightkurve as lk

logging.basicConfig(level=logging.INFO, stream=sys.stdout)
log = logging.getLogger("TESSQuats")
log.setLevel(logging.INFO)

TESSVectors_Products_Base = 'Products'

TESSVectos_local=False
TESSVectors_Local_ENG_Path = 'Eng'
TESSVectors_Local_tpf_fast_path = 'SourceData/fast'
TESSVectors_Local_tpf_short_path = 'SourceData/short/'
TESSVectors_Local_tpf_ffi_path = 'SourceData/FFI/'

def get_ccd_centers(Sector, Camera, CCD):
    """Given a sector/camera/ccd - return the center ra, dec for that TESS pointing
    Uses a subprocess TESSPoint reverse search - could change this to an import?
    """
    # package tess_stars2px.py with this for now
    point = subprocess.run(
        [
            "python",
            "/Users/tapritc2/tessgi/tesspoint/tess-point/tess_stars2px.py",
            "-r",
            str(Sector),
            str(Camera),
            str(CCD),
            str(1024),
            str(1024),
        ],
        capture_output=True,
        text=True,
    )
    ra = float(point.stdout.split(" ")[0])
    dec = float(point.stdout.split(" ")[1])
    return ra, dec


# Get a timing array

def get_times_from_mast(ra, dec, Sector, exptime):
    # This queries external data and occasionally times out
    # Built in some redundancy
    max_attempts = 5
    while max_attempts:
        try:
            if exptime != "ffi":
                # Use TESSCut for an FFI
                results = lk.search_targetpixelfile(
                    f"{ra} {dec}",
                    radius=1e4,
                    sector=Sector,
                    mission="TESS",
                    exptime=exptime,
                    limit=1,
                    author="SPOC",
                ).download(quality_bitmask="none")
            else:
                results = lk.search_tesscut(f"{ra} {dec}", sector=Sector).download(
                    cutout_size=(1, 1), quality_bitmask="none"
                )
            break
        except TimeoutError:
            max_attempts -= 1
    return results

def get_times_local(Sector, Camera, exptime):
    if(exptime != 'ffi'):
        regex = f"*s{Sector:04d}*ffic.fits "
        ffi_list = glob.glob(f"{TESSVectors_Local_tpf_ffi_path}{regex}")
        timing_benchmark = []
        timing_corr = []
        cadences = []
        quality = []
        
        for ffi in ffi_list:
            with fits.open(ffi) as hdu:
                timing_benchmark.append(0.5 * abs(hdu[0].header['TIMESTOP'] -'TIMESTART'))
                timing_corr.append(hdu[1].header['BARYCORR'])
                cadences.append(hdu[0].header['FFIINDEX'])
                quality.append(hdu[1].header['DQUALITY'])
        results = [timing_benchmark, 
                   timing_corr,
                   cadences,
                   quality]
    else:
        if exptime == 'short':
            regex = f"*s{Sector:04d}*s_tp.fits"
            path = TESSVectors_Local_tpf_short_path
        if exptime == 'fast':
            regex = f"*s{Sector:04d}*fast-tp.fits"
            path = TESSVectors_Local_tpf_fast_path
        
        tpf_list = glob.glob(f"{path}{regex}")

        cam = 0
        i=0
        while (cam != Camera) and i < len(tpf_list):
            tpf = lk.read(tpf_list[i])
            tpf.meta['CAMERA']
            i+=1
        if i == len(tpf_list):
            raise RuntimeWarning(f"No tpf for Sector {Sector} Camera {Camera} {exptime} found")
        else:
            results = tpf

    return results

def get_timing_midpoints_tpf(
    ra_c, dec_c, Sector, Camera, exptime, break_times=None, local=False):
    """Given A RA, DEC position, Sector, and observing cadence -
    return an array with the ob served times extracted from a TPF

    This has grown - now also returns FFI filenames and segment breaks as these
    also use information derived from the tpfs
    """

    timing_benchmark = []
    if(exptime != 'ffi' and not (TESSVectos_local or local)):
        # If we are NOT doing FFI Local Processing
        for ra, dec in zip(ra_c, dec_c):
            if (not TESSVectos_local) or (not local):
                results = get_times_from_mast(
                    ra, dec, Sector, exptime
                )
            else:
                results = get_times_local(Sector, Camera, exptime)

            # Get timecorr, correct for it
            # lightkurve PR?
            nonzero_mask = np.nonzero(results.time.value)  

            with fits.open(results.path) as hdu:
                timing_corr = hdu[1].data["TIMECORR"][nonzero_mask]
                timing_ccd = results.time.value[nonzero_mask] - timing_corr
                ffi_list = []
                if exptime == "ffi":
                    ffi_list = hdu[1].data["FFI_FILE"][nonzero_mask]


            if len(timing_benchmark) == 0:
                timing_benchmark = timing_ccd
                cadences = results.cadenceno[nonzero_mask]  # .value
                quality = results.quality[nonzero_mask]
            else:
                if not np.allclose(timing_benchmark, timing_ccd, rtol=1e-6):
                    warnings.warn(f"Timing Cadence Mismatch in {exptime} TPF Data")
                    log.warning(f"Length of Benchmark: {len(timing_benchmark)}")
                    log.warning(f"Length of Comparison CCD: {len(timing_ccd)}")
                    if len(timing_benchmark) == len(timing_ccd):
                        log.warning(
                            f"Maximum Difference between the arrays: {max(abs(timing_benchmark-timing_ccd))}"
                        )
                    # maybe set timing_benchmark = union of two timing sets here
                    # swap above prints to logging-debug
    else:
        # If we are doing ffi local processing
        timing_benchmark, timing_corr, cadences, quality = get_times_local(Sector, Camera, exptime)

    # Calculate break times from 120s data and asign labeling from them
    # if ffi or 020 
    if (break_times is None):
        break_times = get_break_times(timing_benchmark.value)
        
    segment_list = get_segment_label(timing_benchmark.value, break_times)


    # Add in supplementrary Information
    TimingArr = np.zeros((5, len(timing_benchmark)), dtype=np.double)
    TimingArr[0] = timing_benchmark
    TimingArr[1] = timing_corr
    TimingArr[2] = cadences
    TimingArr[3] = quality
    TimingArr[4] = segment_list

    return TimingArr, break_times, ffi_list


def get_camera_sector_cadences(Sector, Camera):
    """For a given Sector, Camera Get the TimingMidpoint Arrays and segmentLabels
    For each Cadence in addition to the FFI File List

    Returns midpoints, segmentslabels, and ffilist arrays
    """
    cutoff_20s = 27  # No 20s data before cycle_whatever
    CCD_List = [1]  # (1,2,3,4)
    # write now we're checking that each CCD has the same timearr for each product
    # once we're done with that - we can can use only CCD 1 or whatever
    # IN PROGRESS
    """ there looks like their can be a small offset on the timing array between CCD's - 
        exposure time? ~1s, using Camera 1 timing as reference [TO DOCUMENT/CONFIRM]"""
    ra_c = []
    dec_c = []
    # Get CCD Centers:
    for CCD in CCD_List:
        ra, dec = get_ccd_centers(Sector, Camera, CCD)
        ra_c.append(ra)
        dec_c.append(dec)

    midpoint_20s = None
    
    # Get 120s TPF Timing
    log.info(f"\t\tGetting 120s Cadence Midpoints Sector: {Sector} Camera: {Camera}")
    midpoint_120s, break_times, _ = get_timing_midpoints_tpf(
        ra_c, dec_c, Sector, Camera, "short"
    )
    
    # update variable name, midpoint is now a list of objects including the midpoints
    if Sector >= cutoff_20s:
        # Get high-cadence 20s timing
        log.info(f"\t\tGetting 20s Cadence Midpoints Sector: {Sector} Camera: {Camera}")
        midpoint_20s, _, _ = get_timing_midpoints_tpf(
            ra_c, dec_c, Sector, Camera, "fast",
            break_times=break_times,
        )

    # Get FFI Timing
    log.info(f"\t\tGetting FFI Cadence Midpoints Sector: {Sector} Camera: {Camera}")
    midpoint_ffi, _, ffi_list = get_timing_midpoints_tpf(
        ra_c,
        dec_c,
        Sector,
        Camera,
        "ffi",
        break_times=break_times,
    )

    midpoints = [midpoint_20s, midpoint_120s, midpoint_ffi]
    return midpoints, ffi_list


def listFD(url, ext=""):
    """gets a list of files from a given URL using the BS4 parcer"""
    page = requests.get(url).text
    soup = BeautifulSoup(page, "html.parser")
    return [
        url + "/" + node.get("href")
        for node in soup.find_all("a")
        if node.get("href").endswith(ext)
    ]

def get_eng_data_local(Sector, eng):
    EngData = glob.glob(f"{TESSVectors_Local_ENG_Path}/*sector{Sector:02d}-{eng}.fits")
    return EngData[0]

def get_eng_data(Sector, eng, local=False):
    log.info(f"\tGetting Engineering Data Type: {eng} for Sector: {Sector}")
    if(not local):
         base_url = "https://archive.stsci.edu/missions/tess/engineering"
         filelist = [file for file in listFD(base_url, f"sector{Sector:02d}-{eng}.fits")]
         # TODO Later should potentially code some logic to use the largest number/latest processing
         EngData = fits.open(filelist[0])
    else:
        EngData = get_eng_data_local(Sector, eng)
    return EngData


def Bin_Quat_Cadence(cadence, CameraData, Camera):
    """Bin engineering quaternions to observed cadences
    cadence - observed cadences / array of observation midpoints
    CameraData - Quaternion Data For Camera from the quat engineering file
    Camera - Camera Number
    """
    ExpTime = np.double(np.median(np.abs(np.diff(cadence))))
    log.debug(f"ExpTime={ExpTime}")
    SubArr_Len = int(ExpTime * (60.0 * 60.0 * 24.0)) # convert days to seconds, 
    # Nbins = ExpTime(seconds), since quats have 2s exposures - this means
    # SubArr_Len is twice as many bins as we expect to need 

    Time_Len = len(cadence)
    Source_Len = len(CameraData["Time"])
    BinnedQuats = np.zeros(
        (17, Time_Len), dtype=np.double
    )  # ( 8 Quat Vectors * 3 + 2 (Tmin/max) + #GS + MSTot + FoM)
    min_ind = int(0)
    max_ind = SubArr_Len  # Days -> Seconds / 2 S per Quat Exp
   
    Flag_ZeroMask = 0
    #Debug Output
    log.debug(f"Index Time Length: {Time_Len}")
    log.debug(f"Source Time Length: {len(CameraData)}")

    for i in range(Time_Len - 1):
        log.debug(f"min: {min_ind} max: {max_ind} max_len: {Time_Len} SubArr_Len: {SubArr_Len}")
        log.debug(f"i: {i} ")

        if i < (Time_Len - 1):
            while CameraData["Time"][max_ind] < cadence[i + 1]:
                max_ind = max_ind + SubArr_Len
        if i > 1:
            while CameraData["Time"][min_ind] > cadence[i - 1]:
                min_ind = min_ind - 1

        SubArr = CameraData[min_ind:max_ind]

        mask = np.nonzero(
            abs(np.double(SubArr["Time"] - cadence[i]))
            < np.double(np.double(0.5) * ExpTime)
        )

        if len(mask[0]) == 0:
            Flag_ZeroMask = Flag_ZeroMask + 1

        #Debug Output
        log.debug(f"\t\t\t\tWarning: 0 Quaternions found in the time bin!")
        log.debug(f"Min: {CameraData['Time'][min_ind]}")
        log.debug(f"Max: {CameraData['Time'][max_ind]}")
        log.debug(f"Cadence: {cadence[i]}")
        log.debug(f"ExpTime: {ExpTime}")
        log.debug(f"ArrEq: {abs(np.double(SubArr['Time'] - cadence[i])) < np.double(np.double(0.5) * ExpTime)}")
        log.debug(f"DMin: {cadence[i] - 0.5 * ExpTime}")
        log.debug(f"DMax: {cadence[i] + 0.5 * ExpTime}")
        log.debug(f"SubArr: {SubArr['Time']}")
        log.debug(f"Abs(SubArr-Cadence): {abs(np.double(SubArr['Time'] - cadence[i]))}")
        log.debug(f"Mask: {mask}")
        log.debug(f"Mask_Len: {len(mask[0])}")
        if len(mask[0]) != 0:
            BinnedQuats[0][i] = np.min(SubArr["Time"][mask])  # Min Midpoint Binned
            BinnedQuats[1][i] = np.max(SubArr["Time"][mask])  # Max Midpoint Binned
            BinnedQuats[2][i] = np.min(
                SubArr[f"C{Camera:d}_FOM"][mask]
            )  # Figure of Merit
            BinnedQuats[3][i] = np.min(
                SubArr[f"C{Camera:d}_NUM_GSUSED"][mask]
            )  # Number of Guide Stars Used # remove this or FoM
            BinnedQuats[4][i] = len(
                mask[0]
            )  # N quats used - diagnose wierdness post binning

            BinnedQuats[5][i] = np.median(SubArr[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[6][i] = np.std(SubArr[f"C{Camera:d}_Q1"][mask])
            BinnedQuats[7][i] = sigma_clipped_stats(SubArr[f"C{Camera:d}_Q1"][mask])[2]

            BinnedQuats[8][i] = np.median(SubArr[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[9][i] = np.std(SubArr[f"C{Camera:d}_Q2"][mask])
            BinnedQuats[10][i] = sigma_clipped_stats(SubArr[f"C{Camera:d}_Q2"][mask])[2]

            BinnedQuats[11][i] = np.median(SubArr[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[12][i] = np.std(SubArr[f"C{Camera:d}_Q3"][mask])
            BinnedQuats[13][i] = sigma_clipped_stats(SubArr[f"C{Camera:d}_Q3"][mask])[2]

            BinnedQuats[14][i] = np.median(SubArr[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[15][i] = np.std(SubArr[f"C{Camera:d}_Q4"][mask])
            BinnedQuats[16][i] = sigma_clipped_stats(SubArr[f"C{Camera:d}_Q4"][mask])[2]

            min_ind = min_ind + int(np.max(mask[0]) + 1)

        else:
            for j in range(0, 16):
                BinnedQuats[j][i] = np.NaN
                # Sector 43 - Camera 2 not being used, only 4
                # Running into issues where there are no quaternions
                # Substitute Nan's?
                # TODO Include data quality flags, check to see if this is in known bad data and
                # throw a harder error?
        # Add Cadence #
        # Quality Masks
        # print(max(mask), SubArr_Len)
        max_ind = int(min_ind) + int(SubArr_Len)
        if max_ind > Source_Len:
            max_ind = Source_Len - 1

    if Flag_ZeroMask > 0:
        log.info(
            f"\t\t\t\tWARNING: Zero Quaternions found for {Flag_ZeroMask} observed cadences"
        )
    return BinnedQuats


def EMI_Extension_Dict(ext_string):
    """Dictionary of fits emi file table column names"""
    # This could be done algorithmically, but thought it was less error prone this way
    # I broke this out of the function b/c it was big & ugly
    emi_dict = {
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
        "emi.cam4_moon_azimuth": 22,
    }
    return emi_dict[ext_string]


def Bin_EMI_Cadence(cadence, EMIData, Camera, type):
    """Bin the EMI information to the observed cadences
    cadence - observed time midpoint array / end-product observation cadences
    EMIData - EMI data file from the emi engineering file
    Camera - camera number
    type - observing cadence for cadence
    """
    typedict = {1: "020", 2: "120", 3: "FFI"}
    ExpTime = np.double(
        abs(cadence[1] - cadence[0])
    )  # Should this be hard-coded instead? prbly med
    if isinstance(type, int):
        type = typedict[type]
    match type:
        case "020":
            SubArr_Len = int(4)
        case "120":
            SubArr_Len = int(4)
        case "FFI":
            SubArr_Len = int(
                (ExpTime / (EMIData[1].data["TIME"][1] - EMIData[1].data["TIME"][0]))
                * 4
            )
        case _:
            log.info("Error - Bin_EMI_Cadence can't match cadence type")
            # throw a logging warning instead

    Time_Len = len(cadence)
    Source_Len = len(EMIData[1].data["Time"])
    BinnedEMI = np.zeros(
        (6, Time_Len), dtype=np.double
    )  # ( (Earth + Moon) ** (Distance + Elev. + Azim.)
    min_ind = int(0)
    max_ind = int(SubArr_Len)

    for i in range(Time_Len - 1):
        if i < (Time_Len - 2):
            while EMIData[1].data["Time"][max_ind] < cadence[i + 1]:
                max_ind = max_ind + SubArr_Len

        if i >= 1:
            while EMIData[1].data["Time"][min_ind] > cadence[i - 1]:
                min_ind = min_ind - 1

        SubArr = EMIData[1].data["Time"][min_ind:max_ind]
        match_ind = np.argmin(abs(np.double(SubArr - cadence[i])))
        # print(f"i: {i} min_ind: {min_ind} max_ind:{max_ind} match_ind:{match_ind} Time_Len:{Time_Len} Subarr_Len:{SubArr_Len} ")

        if isinstance(min_ind, int):
            # Earth Distance
            BinnedEMI[0][i] = EMIData[EMI_Extension_Dict("emi.earth_distance")].data[
                min_ind:max_ind
            ]["VALUE"][match_ind]
            # Earth Camera Elev
            BinnedEMI[1][i] = EMIData[
                EMI_Extension_Dict(f"emi.cam{Camera}_earth_boresight_angle")
            ].data[min_ind:max_ind]["VALUE"][match_ind]
            # EarthCamera Azim
            BinnedEMI[2][i] = EMIData[
                EMI_Extension_Dict(f"emi.cam{Camera}_earth_azimuth")
            ].data[min_ind:max_ind]["VALUE"][match_ind]
            # Moon Distance
            BinnedEMI[3][i] = EMIData[EMI_Extension_Dict("emi.moon_distance")].data[
                min_ind:max_ind
            ]["VALUE"][match_ind]
            # Moon Camera Elev
            BinnedEMI[4][i] = EMIData[
                EMI_Extension_Dict(f"emi.cam{Camera}_moon_boresight_angle")
            ].data[min_ind:max_ind]["VALUE"][match_ind]
            # Moon Camera Azim
            BinnedEMI[5][i] = EMIData[
                EMI_Extension_Dict(f"emi.cam{Camera}_moon_azimuth")
            ].data[min_ind:max_ind]["VALUE"][match_ind]
        else:
            log.info("Bin_EMI_Cadence can't match a datapoint")
            # throw a real logging.warning

        min_ind = min_ind + int(match_ind)
        max_ind = int(min_ind) + int(SubArr_Len)
        if max_ind > Source_Len:
            max_ind = Source_Len - 1

    return BinnedEMI


def Bin_Quat_Camera(CameraData, TimeArr, Camera):
    """Bin quaternions for all cadences for a given camera
    CameraData - quaternion data from the quat engineering file
    TimeArr - observed cadence midpoints from tpfs
    Camera - Camera Number
    """
    Bin20 = None
    if TimeArr[0] is not None:
        log.info("\t\tBinning 20s Quaternion Data")
        Bin20 = Bin_Quat_Cadence(TimeArr[0][0], CameraData, Camera)
    log.info("\t\tBinning 120s Quaternion Data")
    Bin120 = Bin_Quat_Cadence(TimeArr[1][0], CameraData, Camera)
    log.info("\t\tBinning FFI Quaternion Data")
    BinFFI = Bin_Quat_Cadence(TimeArr[2][0], CameraData, Camera)

    return Bin20, Bin120, BinFFI


def Bin_EMI_Camera(EMIData, TimeArr, Camera):
    """Bin emi info for all cadences for a given camera
    EMIData - EMI data from the emi engineering file
    TimeArr - observed cadence midpoints from tpfs
    Camera - Camera Number
    """
    EMI20 = None
    if TimeArr[0] is not None:
        log.info("\t\tBinning 20s EMI Data")
        EMI20 = Bin_EMI_Cadence(TimeArr[0][0], EMIData, Camera, "020")
    log.info("\t\tBinning 120s EMI Data")
    EMI120 = Bin_EMI_Cadence(TimeArr[1][0], EMIData, Camera, "120")
    log.info("\t\tBinning FFI EMI Data")
    EMIFFI = Bin_EMI_Cadence(TimeArr[2][0], EMIData, Camera, "FFI")

    return EMI20, EMI120, EMIFFI

def get_break_times(times):
    break_inds = np.where(np.diff(np.asarray(times)) > 1/24)[0] + 1
    break_times = np.asarray([np.mean(times[ind -1:ind+1]) for ind in break_inds])
    return break_times

def get_segment_label(times, break_times):
    """Create labels for segments - using the quality flag and timing array

    right now, an 'segment' is every stretch of observations between earth point
    quality flags (=8), with each earth point stretch causing an increment in
    segment label
    """
    label = 1
    # cut should not be nescessary
    cut_times=times[times != 0]
    segment_labels = np.zeros_like(cut_times)

    for btime in break_times:
        if(label == 1):
            segment_labels[cut_times < btime] = label
        label += 1
        segment_labels[cut_times > btime] = label
    if(np.any(segment_labels == 0)):
       raise RuntimeWarning("get_segment_label warning: Not all times have a segment label")

    return segment_labels

def write_vector_sector_camera(
    BinnedQuats, BinnedEMI, FFIList, TimeArr, Sector, Camera
):
    """Write TESSVectors Info to CSV files for all observed cadences"""
    typedict = {1: "020", 2: "120", 3: "FFI"}
    #Binned_Dir = "TESSVectors_products"
    for Quat, EMI, Time, i in zip(
        BinnedQuats, BinnedEMI, TimeArr, [1, 2, 3]
    ):
        from datetime import date

        fname = f"{TESSVectors_Products_Base}/Vectors/{typedict[i]}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{typedict[i]}.csv"
        log.info(
            f"\t\tWriting Sector: {Sector} Camera: {Camera} Cadence {typedict[i]} to:"
        )
        log.info(f"\t\t\t to {fname}")

        if Time is not None:
            # First write the File Header
            # this will over-write by default
            with open(fname, "w") as f:
                f.write(f"# TESS Quaternions downsampled to end-user cadences\n")
                f.write(f"# Sector: \t{Sector}\n")
                f.write(f"# Camera: \t{Camera}\n")
                f.write(f"# Cadence:\t{typedict[i]}\n")
                f.write(
                    f"# This file contains TESS quaternions that have been downsampled from the TESS\n"
                )
                f.write(
                    f"# native 2-second exposures to the end-user product cadence (e.g .20s/120s/FFI)\n"
                )
                f.write(f"# For more information See:\n")
                f.write(
                    f"#     - The github repo that created this file at https://github.com/tylerapritchard/TESSQuats\n"
                )  # This will presumably move to tessgi at some point
                f.write(
                    f"#     - The TESS Instrument Handbook at https://archive.stsci.edu/missions-and-data/tess\n"
                )
                f.write(
                    f"#     - The original quaternion engineering files at https://archive.stsci.edu/missions-and-data/tess/data-products\n"
                )
                f.write(
                    f"# Please check the TESS Sector Data Release Notes at https://archive.stsci.edu/tess/tess_drn.html\n"
                )
                f.write(
                    f"# for information regarding the telescope pointing and tracking (in addition to other image quality information) as\n"
                )
                f.write(
                    f"# for some sectors/cameras alternative camera quaternions are substituted for a variety of programattic reasons.\n\n"
                )
                f.write(f"# Column Descriptions:\n")
                f.write(f"# Cadence #: Cadence index from the source tpf\n")
                f.write(
                    f"# MidTime: The exposure midpoint in spacecraft corrected time (i.e. tpf.time)\n"
                )
                f.write(
                    f"# TimeCorr: The Time Correction for spacecraft to Barycentric time at that cadence\n"
                )
                f.write(
                    f"# Quat_Start: The timestamp of the earliest quaternion used in the bin\n"
                )
                f.write(
                    f"# Quat_Stop: The timestamp of the last quaternion used in the bin\n"
                )
                f.write(
                    f"# Quat_MIN_FOM: The worst Figure of Merit from the source quaternions\n"
                )
                f.write(
                    f"# Quat_MIN_NUM_GSUSED: The lowest number of guide stars used in the source quaternions\n"
                )
                f.write(
                    f"# Quat_NBinned: The number of quaternions binned into this final result.\n"
                )
                f.write(
                    f"# Quat[1-4]_Med: The Quaternion #[1-4] median value from the binned values \n"
                )
                f.write(
                    f"# Quat[1-4]_StdDev: The standard deviation of Quaternion #[1-4] binned values\n"
                )
                f.write(
                    f"# Quat[1-4]_SigClip: The Sigma-Clipped Standard Deviation of Quaternion #[1-4] binned values\n"
                )
                f.write(
                    f"# Quat[1-4]_CRM_Med: The Quaternion #[1-4] median value with the highest and lowest values excluded \n\n"
                )
                f.write(f"# Earth_Distance: Distance to Earth in Earth Radii \n\n")
                f.write(
                    f"# Earth_Camera_Angle: Angle of Earth from Camera Boresight in Degrees \n\n"
                )
                f.write(
                    f"# Earth_Camera_Azimuth: Azimuth of Earth around Camera Boresight in Degrees \n\n"
                )
                f.write(f"# Moon_Distance: Distance to Moon in Earth Radii \n\n")
                f.write(
                    f"# Moon_Camera_Angle: Angle of Moon from Camera Boresight in Degrees \n\n"
                )
                f.write(
                    f"# Moon_Camera_Azimuth: Azimuth of Moon around Camera Boresight in Degrees \n\n"
                )
                if i == 3:
                    f.write(
                        f"# FFIFile: The FFI file assosciated with this observing cadence\n"
                    )
                f.write(f"# Processing Date-Time: {date.today()}\n\n")
            df = pd.DataFrame(
                data={
                    "Cadence": Time[2],
                    "MidTime": Time[0] + Time[1],
                    "TimeCorr": Time[1],
                    "Segment": Time[4],
                    "Quat_Start": Quat[0],
                    "Quat_Stop": Quat[1],
                    "Quat_MIN_FOM": Quat[2],
                    "Quat_MIN_NUM_GSUSED": Quat[3],
                    "Quat_NBinned": Quat[4],
                    "Quat1_Med": Quat[5],
                    "Quat1_StdDev": Quat[6],
                    "Quat1_StdDev_SigClip": Quat[7],
                    "Quat2_Med": Quat[8],
                    "Quat2_StdDev": Quat[9],
                    "Quat2_StdDev_SigClip": Quat[10],
                    "Quat3_Med": Quat[11],
                    "Quat3_StdDev": Quat[12],
                    "Quat3_StdDev_SigClip": Quat[13],
                    "Quat4_Med": Quat[14],
                    "Quat4_StdDev": Quat[15],
                    "Quat4_StdDev_SigClip": Quat[16],
                    "Earth_Distance": EMI[0],
                    "Earth_Camera_Angle": EMI[1],
                    "Earth_Camera_Azimuth": EMI[2],
                    "Moon_Distance": EMI[3],
                    "Moon_Camera_Angle": EMI[4],
                    "Moon_Camera_Azimuth": EMI[5],
                }
            )
            if i == 3:
                df["FFIFile"] = FFIList

            # df.astype({'Cadence':int,
            #           'Sector': int,
            #           'Camera':int}
            #           #'Quat_MIN_NUM_GSUSED':int,
            #           #'Quat_NBinned': int}
            #        ).to_csv(fname, index=False, mode = "a")
            df.to_csv(fname, index=False, mode="a")

def create_vectors_sector(Sector):
    """For a given sector, create TESSVectors Information and write to a CSV file"""

    log.info(f"Starting Sector: {Sector}")
    QuatData = get_eng_data(Sector, 'quat')
    EMIData = get_eng_data(Sector, 'emi')

    for Camera in [1, 2, 3, 4]:
        log.info(f"\tStarting Camera: {Camera}")
        TimeArr, FFIList = get_camera_sector_cadences(Sector, Camera)
        BinnedQuats = Bin_Quat_Camera(QuatData[Camera].data, TimeArr, Camera)
        BinnedEMI = Bin_EMI_Camera(EMIData, TimeArr, Camera)
        write_vector_sector_camera(
            BinnedQuats, BinnedEMI, FFIList, TimeArr, Sector, Camera
        )


def quality_to_color(qual):
    default = TessQualityFlags.create_quality_mask(qual, bitmask="default")
    hard = TessQualityFlags.create_quality_mask(qual, bitmask="hard")
    # hardest = TessQualityFlags.create_quality_mask(qual, bitmask='hardest')

    carr = np.zeros(len(qual))
    # carr[~ hardest] = 0.33
    carr[~hard] = 0.5
    carr[~default] = 1

    return carr


def plot_quat(axs, time, quat, dev, qual, QuatLabel):
    nsigma = 1

    # norm = np.nanmedian(np.double(quat[~np.isfinite(quat)]))
    norm_quat = quat / norm
    norm = np.nanmedian(quat)
    norm_dev = dev / norm
    # mean_norm_dev = np.nanmedian(dev[~np.isfinite(dev)] / norm)
    mean_norm_dev = np.nanmedian(dev / norm)
    tmin = np.min(time)
    tmax = np.max(time)

    ymin = np.nanmedian(norm_quat) - 3 * mean_norm_dev
    ymax = np.nanmedian(norm_quat) + 3 * mean_norm_dev

    im = axs.imshow(
        np.vstack((qual,)),
        extent=(tmin, tmax, ymin, ymax),
        interpolation="nearest",
        aspect="auto",
        cmap=cm.PuRd,
        vmax=1,
    )

    axs.scatter(time, norm_quat, s=3, label="Median Quaternion Value", color="k")

    axs.fill_between(
        time,
        norm_quat - (nsigma * norm_dev),
        norm_quat + (nsigma * norm_dev),
        alpha=0.6,
        color="grey",
    )

    axs.set_xlim(tmin, tmax)
    axs.set_ylim(ymin, ymax)

    axs.tick_params(axis="x", labelsize=0)

    axs.set_ylabel(f"{QuatLabel}", weight="bold", size=24)
    axs.set_yticks([])

    return im


def create_diagnostic_timeseries(Sector, Camera, Cadence):
    # TODO add segment labels to plot?
    typedict = {1: "020", 2: "120", 3: "FFI"}
    if type(Cadence) != str:
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    fname = f"{TESSVectors_Products_Base}/Vectors/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"

    nplots = 3
    if os.path.isfile(fname):
        quatdf = pd.read_csv(fname, comment="#", index_col=False)
        qual_im = quality_to_color(quatdf.Quality)

        fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))
        im = plot_quat(
            axs[0],
            quatdf.MidTime,
            quatdf.Quat1_Med,
            quatdf.Quat1_StdDev,
            qual_im,
            "Quaternion 1",
        )
        im = plot_quat(
            axs[1],
            quatdf.MidTime,
            quatdf.Quat2_Med,
            quatdf.Quat2_StdDev,
            qual_im,
            "Quaternion 2",
        )
        im = plot_quat(
            axs[2],
            quatdf.MidTime,
            quatdf.Quat3_Med,
            quatdf.Quat3_StdDev,
            qual_im,
            "Quaternion 3",
        )
        plt.subplots_adjust(hspace=0)
        axs[0].set_title(
            f"TESS Sector {Sector} Camera {Camera} Quaternions", weight="bold", size=26
        )
        axs[-1].tick_params(axis="x", labelsize=18)
        axs[-1].set_xlabel("TESS BTJD", weight="bold", size=24)

        cax = plt.axes([0.92, 0.11, 0.075, 0.77])
        cbar = plt.colorbar(mappable=im, cax=cax, ticks=[0, 0.5, 1])
        cbar.ax.set_yticklabels(["Unflagged", "Aggressive", "Conservative"], size=18)
        cbar.set_label("Data Flagging Level (Lower Is Better)", size=24, weight="bold")

        # plt.tight_layout()
        fout = f"{TESSVectors_Products_Base}/Diagnostics/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_Quat.png"
        plt.show()
        plt.savefig(fout, dpi=300, bbox_inches="tight", rasterize=True)
        plt.close(fig)


def plot_lsperiodogram(ax, time, median, std, QuatLabel):
    lc = lk.LightCurve(data={"time": time, "flux": std})
    ls = lc.to_periodogram(maximum_period=13.5)
    ls.plot(ax=ax, lw=0.1, color="k", ylabel=" ")

    ax.set_ylabel(f"{QuatLabel} Power", weight="bold", size=24)
    # ax.set_yticks([])
    ax.tick_params(axis="y", labelsize=18)
    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_ylim(min(ls.power), max(ls.power))

    return


def create_diagnostic_periodogram(Sector, Camera):
    Cadence = 1
    typedict = typedict = {1: "020", 2: "120", 3: "FFI"}
    if type(Cadence) != str:
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    Binned_Dir = "TESSVectors_products"
    fname = f"{TESSVectors_Products_Base}/Vectors/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"

    nplots = 3
    if os.path.isfile(fname):
        quatdf = pd.read_csv(fname, comment="#", index_col=False)
        fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))
        plot_lsperiodogram(
            axs[0],
            quatdf.MidTime,
            quatdf.Quat1_Med,
            quatdf.Quat1_StdDev,
            "Quaternion 1",
        )
        plot_lsperiodogram(
            axs[1],
            quatdf.MidTime,
            quatdf.Quat2_Med,
            quatdf.Quat2_StdDev,
            "Quaternion 2",
        )
        plot_lsperiodogram(
            axs[2],
            quatdf.MidTime,
            quatdf.Quat3_Med,
            quatdf.Quat3_StdDev,
            "Quaternion 3",
        )
        plt.subplots_adjust(hspace=0)

        axs[0].set_title(
            f"TESS Sector {Sector} Camera {Camera} Quaternion Power Spectra",
            weight="bold",
            size=26,
        )
        axs[-1].tick_params(axis="x", labelsize=18)
        axs[-1].set_xlabel("Period [days]", weight="bold", size=24)

        forward = lambda x: x * 24.0 * 3600.0
        inverse = lambda x: x / 24.0 / 3600.0
        ax2 = axs[0].secondary_xaxis("top", functions=(forward, inverse))
        ax2.set_xlabel("Period [seconds]", weight="bold", size=24)
        ax2.tick_params(axis="x", labelsize=18)

        fout = f"{TESSVectors_Products_Base}/Diagnostics/Periodograms/TessVectors_S{Sector:03d}_C{Camera}_QuatPower.png"
        plt.savefig(fout, dpi=300, bbox_inches="tight")
        plt.close(fig)


def create_diagnostic_emi(Sector, Camera, Cadence):
    import matplotlib as mpl

    mpl.rcParams["agg.path.chunksize"] = 10000000

    typedict = typedict = {1: "020", 2: "120", 3: "FFI"}
    if type(Cadence) != str:
        cadence_name = typedict[Cadence]
    else:
        cadence_name = Cadence
    fname = f"{TESSVectors_Products_Base}/Vectors/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}.csv"

    nplots = 3
    if os.path.isfile(fname):
        quatdf = pd.read_csv(fname, comment="#", index_col=False)

        qual_im = quality_to_color(quatdf.Quality)
        fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))

        axs[0].scatter(
            quatdf.MidTime, quatdf.Earth_Distance, label="Earth", color="seagreen"
        )
        axs[0].scatter(
            quatdf.MidTime, quatdf.Moon_Distance, label="Moon", color="darkturquoise"
        )
        axs[0].legend(
            prop={"weight": "bold", "size": 24}, scatterpoints=3, markerscale=2
        )

        tmin = np.min(quatdf.MidTime)
        tmax = np.max(quatdf.MidTime)

        ymin = min(min(quatdf.Earth_Distance), min(quatdf.Moon_Distance))
        ymax = max(max(quatdf.Earth_Distance), max(quatdf.Moon_Distance))
        im = axs[0].imshow(
            np.vstack((qual_im,)),
            extent=(tmin, tmax, ymin, ymax),
            interpolation="nearest",
            aspect="auto",
            cmap=cm.PuRd,
            vmax=1,
        )

        axs[1].scatter(
            quatdf.MidTime, quatdf.Earth_Camera_Angle, label="Earth", color="seagreen"
        )
        axs[1].scatter(
            quatdf.MidTime,
            quatdf.Moon_Camera_Angle,
            label="Moon",
            color="darkturquoise",
        )
        axs[1].legend(
            prop={"weight": "bold", "size": 24}, scatterpoints=3, markerscale=2
        )

        ymin = min(min(quatdf.Earth_Camera_Angle), min(quatdf.Moon_Camera_Angle))
        ymax = max(max(quatdf.Earth_Camera_Angle), max(quatdf.Moon_Camera_Angle))
        im = axs[1].imshow(
            np.vstack((qual_im,)),
            extent=(tmin, tmax, ymin, ymax),
            interpolation="nearest",
            aspect="auto",
            cmap=cm.PuRd,
            vmax=1,
        )

        axs[2].scatter(
            quatdf.MidTime, quatdf.Earth_Camera_Azimuth, label="Earth", color="seagreen"
        )
        axs[2].scatter(
            quatdf.MidTime,
            quatdf.Moon_Camera_Azimuth,
            label="Moon",
            color="darkturquoise",
        )
        axs[2].legend(
            prop={"weight": "bold", "size": 24}, scatterpoints=3, markerscale=2
        )

        ymin = min(min(quatdf.Earth_Camera_Azimuth), min(quatdf.Moon_Camera_Azimuth))
        ymax = max(max(quatdf.Earth_Camera_Azimuth), max(quatdf.Moon_Camera_Azimuth))
        im = axs[2].imshow(
            np.vstack((qual_im,)),
            extent=(tmin, tmax, ymin, ymax),
            interpolation="nearest",
            aspect="auto",
            cmap=cm.PuRd,
            vmax=1,
        )

        axs[0].set_ylabel("Distance", weight="bold", size=24)
        axs[1].set_ylabel("Camera {Camera} Angle", weight="bold", size=24)
        axs[2].set_ylabel("Camera {Camera} Azimuth", weight="bold", size=24)

        axs[2].set_xlabel("TESS BTJD", weight="bold", size=24)

        axs[0].tick_params(axis="y", labelsize=18)
        axs[1].tick_params(axis="y", labelsize=18)
        axs[2].tick_params(axis="y", labelsize=18)

        axs[2].tick_params(axis="x", labelsize=18)
        plt.subplots_adjust(hspace=0)

        cax = plt.axes([0.92, 0.11, 0.075, 0.77])
        cbar = plt.colorbar(mappable=im, cax=cax, ticks=[0, 0.5, 1])
        cbar.ax.set_yticklabels(["Unflagged", "Aggressive", "Conservative"], size=18)
        cbar.set_label("Data Flagging Level (Lower Is Better)", size=24, weight="bold")

        fout = f"{TESSVectors_Products_Base}/Diagnostics/{cadence_name}_Cadence/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_emi.png"
        plt.savefig(fout, dpi=300, bbox_inches="tight")
        # plt.show()
        plt.close(fig)


def create_diagnostics_sector(Sector):
    # CameraCadence
    # Sector, Camera, Cadence = SectorCameraCadence
    log.info(f"Creating Diagnostics for Sector: {Sector}")
    for Camera in [1, 2, 3, 4]:
        create_diagnostic_periodogram(Sector, Camera)
        for Cadence in [1, 2, 3]:
            create_diagnostic_timeseries(Sector, Camera, Cadence)
            # Should I create the periodograms from the "raw" 2s data?  probably?
            create_diagnostic_emi(Sector, Camera, Cadence)



def TESSVectors_process_sector(Sector):
    create_vectors_sector(Sector)
    try:
        create_diagnostics_sector(Sector)
    except:
        log.warning("\t\t\t Warning, Plotting failed")  # add a real warning
        pass

    return (Sector, True)


def run_bulk_diagnostics(sector_min = 1, sector_max = 65, camera_min = 1, camera_max = 4, processes=7):
    sector_list = range(sector_min, sector_max)
    camera_list = range(camera_min, camera_max)
    cadence_list = range(1, 3)
    inlist = list(product(sector_list, camera_list, cadence_list))

    from multiprocessing.pool import Pool
    from itertools import product

    pool = Pool(processes=processes)
    res = []
    for result in pool.map(create_diagnostics_bulk, inlist):
        res = [res, result]
    pool.close()


def run_bulk_vectors(sector_min = 1, sector_max = 69, processes=7):
    if not processes:
        processes = 7
    sector_list = range(sector_min, sector_max)

    from multiprocessing.pool import Pool

    pool = Pool(processes=processes)
    res = []
    for result in pool.map(create_vectors_sector, sector_list):
        res = [res, result]
    pool.close()


def run_bulk_processing(sector_min = 1, sector_max = 65):
    from multiprocessing.pool import Pool
    import multiprocessing as mp

    sector_list = range(sector_min, sector_max)

    pool = Pool(processes=mp.cpu_count())
    res = []
    for result in pool.map(TESSVectors_process_sector, sector_list):
        res = [res, result]
    pool.close()
