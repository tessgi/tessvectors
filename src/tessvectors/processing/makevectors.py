import os
import sys
import glob
import subprocess
import warnings
import logging
import requests

import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

from lightkurve import TessQualityFlags
import numpy as np

from bs4 import BeautifulSoup
import pandas as pd

from itertools import product
from datetime import date

from copy import deepcopy

mpl.rcParams["agg.path.chunksize"] = 10000000

import lightkurve as lk

logging.basicConfig(level=logging.INFO, stream=sys.stdout)
log = logging.getLogger("tessvectors")
log.setLevel(logging.INFO)

# logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
# log = logging.getLogger("TESSQuats")
# log.setLevel(logging.DEBUG)


class makevectors(object):
    def __init__(
        self,
        TESSVectors_Products_Base="Products",
        TESSVectos_local=False,
        TESSVectors_Local_ENG_Path="Eng",
        TESSVectors_Local_tpf_fast_path="SourceData/fast",
        TESSVectors_Local_tpf_short_path="SourceData/short/",
        TESSVectors_Local_tpf_ffi_path="SourceData/FFI/",
        check_exists=True,
    ):
        self.TESSVectors_Products_Base = TESSVectors_Products_Base
        self.TESSVectos_local = TESSVectos_local
        self.TESSVectors_Local_ENG_Path = TESSVectors_Local_ENG_Path
        self.TESSVectors_Local_tpf_fast_path = TESSVectors_Local_tpf_fast_path
        self.TESSVectors_Local_tpf_short_path = TESSVectors_Local_tpf_short_path
        self.TESSVectors_Local_tpf_ffi_path = TESSVectors_Local_tpf_ffi_path
        self.check_exists = check_exists

        self.make_dir_structure()

        self.typedict = {1: "020", 2: "120", 3: "FFI"}

    def make_dir_structure(self):
        os.makedirs(self.TESSVectors_Products_Base, exist_ok=True)
        os.makedirs(self.TESSVectors_Products_Base + "/Vectors", exist_ok=True)
        os.makedirs(
            self.TESSVectors_Products_Base + "/Vectors/020_Cadence", exist_ok=True
        )
        os.makedirs(
            self.TESSVectors_Products_Base + "/Vectors/120_Cadence", exist_ok=True
        )
        os.makedirs(
            self.TESSVectors_Products_Base + "/Vectors/FFI_Cadence", exist_ok=True
        )
        os.makedirs(self.TESSVectors_Products_Base + "/Diagnostics", exist_ok=True)
        os.makedirs(
            self.TESSVectors_Products_Base + "/Diagnostics/Periodograms", exist_ok=True
        )
        os.makedirs(self.TESSVectors_Products_Base + "/Diagnostics/EMI", exist_ok=True)
        os.makedirs(
            self.TESSVectors_Products_Base + "/Diagnostics/Vectors", exist_ok=True
        )

    def _vector_base_file(self, Cadence, Sector, Camera):
        return f"TessVectors_S{Sector:03d}_C{Camera}_{Cadence}"

    def _vector_file(self, Cadence, Sector, Camera):
        base_file = self._vector_base_file(Cadence, Sector, Camera)
        return f"{self.TESSVectors_Products_Base}/Vectors/{Cadence}_Cadence/{base_file}"

    def get_ccd_centers(self, Sector, Camera, CCD):
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

    def get_times_from_mast(self, ra, dec, Sector, exptime):
        # This queries external data and occasionally times out
        # Built in some redundancy
        max_attempts = 5
        while max_attempts:
            try:
                if exptime != "ffi":
                    # Use TESSCut for an FFI
                    # results = lk.search_targetpixelfile(
                    #    f"{ra} {dec}",
                    #    radius=1e4,
                    #    sector=Sector,
                    #    mission="TESS",
                    #    exptime=exptime,
                    #    limit=1,
                    #    author="SPOC",
                    # ).download(quality_bitmask="none")
                    if exptime == "short":
                        script_link = f"https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{Sector}_tp.sh"
                    if exptime == "fast":
                        script_link = f"https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{Sector}_fast-tp.sh"
                    scripts = pd.read_csv(
                        script_link, comment="#", sep=" ", usecols=[6], names=["link"]
                    )

                    cam = -1
                    i = -1
                    while cam != 1:
                        i += 1
                        results = lk.TessTargetPixelFile(
                            scripts["link"].values[i], quality_bitmask="none"
                        )
                        cam = results.hdu[0].header["CCD"]

                else:
                    results = lk.search_tesscut(f"{ra} {dec}", sector=Sector).download(
                        cutout_size=(1, 1), quality_bitmask="none"
                    )
                break
            except TimeoutError:
                max_attempts -= 1
        return results

    def get_times_local(self, Sector, Camera, exptime):
        if exptime == "ffi":
            # We need to choose 1 ccd per camera to index our times off of, here we'll use CCD 1
            ccd_index = 1
            regex = f"*s{Sector:04d}-{Camera}-{ccd_index}*ffic.fits "
            ffi_list = glob.glob(f"{self.TESSVectors_Local_tpf_ffi_path}{regex}")
            timing_benchmark = []
            timing_corr = []
            cadences = []
            quality = []

            for ffi in ffi_list:
                with fits.open(ffi) as hdu:
                    timing_benchmark.append(
                        0.5 * abs(hdu[0].header["TIMESTOP"] - "TIMESTART")
                    )
                    timing_corr.append(hdu[1].header["BARYCORR"])
                    cadences.append(hdu[0].header["FFIINDEX"])
                    quality.append(hdu[1].header["DQUALITY"])
            results = [timing_benchmark, timing_corr, cadences, quality]
        else:
            if exptime == "short":
                regex = f"*s{Sector:04d}*s_tp.fits"
                path = self.TESSVectors_Local_tpf_short_path

            if exptime == "fast":
                regex = f"*s{Sector:04d}*fast-tp.fits"
                path = self.TESSVectors_Local_tpf_fast_path
            tpf_list = glob.glob(f"{path}{regex}")

            cam = 0
            i = 0
            while (cam != Camera) and i < len(tpf_list):
                tpf = lk.read(tpf_list[i])
                tpf.meta["CAMERA"]
                i += 1
            if i == len(tpf_list):
                raise RuntimeWarning(
                    f"No tpf for Sector {Sector} Camera {Camera} {exptime} found"
                )
            else:
                results = tpf

        return results

    def get_timing_midpoints_tpf(
        self, ra_c, dec_c, Sector, Camera, exptime, break_times=None, local=False
    ):
        """Given A RA, DEC position, Sector, and observing cadence -
        return an array with the ob served times extracted from a TPF

        This has grown - now also returns FFI filenames and segment breaks as these
        also use information derived from the tpfs
        """
        local_processing = self.TESSVectos_local or local
        timing_benchmark = []
        if (exptime != "ffi") or ((exptime == "ffi") and (not local_processing)):
            # If we are NOT doing FFI Local Processing
            for ra, dec in zip(ra_c, dec_c):
                if not local_processing:
                    log.info(
                        f"\t\t\tUsing Remote Data To Retrieve Times Sector:{Sector} Camera: {Camera} ExpTime: {exptime}"
                    )
                    results = self.get_times_from_mast(ra, dec, Sector, exptime)
                else:
                    log.info(
                        f"\t\t\tUsing LOCAL Data To Retrieve Times Sector:{Sector} Camera: {Camera} ExpTime: {exptime}"
                    )
                    results = self.get_times_local(Sector, Camera, exptime)

                # Get timecorr, correct for it
                # lightkurve PR?
                nonzero_mask = np.nonzero(results.time.value)

                with fits.open(results.path) as hdu:
                    timing_corr = hdu[1].data["TIMECORR"][nonzero_mask]
                    timing_ccd = results.time.value[nonzero_mask] - timing_corr
                    ffi_list = []
                    if exptime == "ffi":
                        # We're essentially grabbing a random tesscut tpf to use for FFI indexing
                        # Here we swap our random camera with the '{camera}' string and add the extension
                        # for a calibrated ffi
                        ffi_list = hdu[1].data["FFI_FILE"][nonzero_mask]
                        ffi_list = [
                            f[0:-12] + "{ccd}" + f[-11:] + "c.fits" for f in ffi_list
                        ]

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
            log.info(
                f"\t\t\tUsing LOCAL Data To Retrieve Times Sector:{Sector} Camera: {Camera} ExpTime: {exptime}"
            )
            timing_benchmark, timing_corr, cadences, quality = self.get_times_local(
                Sector, Camera, exptime
            )

        # Calculate break times from 120s data and asign labeling from them
        # if ffi or 020
        if break_times is None:
            break_times = self.get_break_times(timing_benchmark)

        segment_list = self.get_segment_label(timing_benchmark, break_times)

        # Add in supplementrary Information
        TimingArr = np.zeros((5, len(timing_benchmark)), dtype=np.double)
        TimingArr[0] = timing_benchmark
        TimingArr[1] = timing_corr
        TimingArr[2] = cadences
        TimingArr[3] = quality
        TimingArr[4] = segment_list

        return TimingArr, break_times, ffi_list

    def get_camera_sector_cadences(self, Sector, Camera):
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
            ra, dec = self.get_ccd_centers(Sector, Camera, CCD)
            ra_c.append(ra)
            dec_c.append(dec)

        midpoint_20s = None

        # Get 120s TPF Timing
        log.info(
            f"\t\tGetting 120s Cadence Midpoints Sector: {Sector} Camera: {Camera}"
        )
        midpoint_120s, break_times, _ = self.get_timing_midpoints_tpf(
            ra_c, dec_c, Sector, Camera, "short"
        )

        # update variable name, midpoint is now a list of objects including the midpoints
        if Sector >= cutoff_20s:
            # Get high-cadence 20s timing
            log.info(
                f"\t\tGetting 20s Cadence Midpoints Sector: {Sector} Camera: {Camera}"
            )
            midpoint_20s, _, _ = self.get_timing_midpoints_tpf(
                ra_c,
                dec_c,
                Sector,
                Camera,
                "fast",
                break_times=break_times,
            )

        # Get FFI Timing
        log.info(f"\t\tGetting FFI Cadence Midpoints Sector: {Sector} Camera: {Camera}")
        midpoint_ffi, _, ffi_list = self.get_timing_midpoints_tpf(
            ra_c,
            dec_c,
            Sector,
            Camera,
            "ffi",
            break_times=break_times,
        )

        midpoints = [midpoint_20s, midpoint_120s, midpoint_ffi]
        return midpoints, ffi_list

    def listFD(self, url, ext=""):
        """gets a list of files from a given URL using the BS4 parcer"""
        page = requests.get(url).text
        soup = BeautifulSoup(page, "html.parser")
        return [
            url + "/" + node.get("href")
            for node in soup.find_all("a")
            if node.get("href").endswith(ext)
        ]

    def get_eng_data_local(Sector, eng):
        # something of a placeholder
        # assumes engineering files are being stored in a bulk directory somewhere
        EngData = glob.glob(
            f"{self.TESSVectors_Local_ENG_Path}/*sector{Sector:02d}-{eng}.fits"
        )
        return EngData[0]

    def get_eng_data(self, Sector, eng, local=False):
        log.info(f"\tGetting Engineering Data Type: {eng} for Sector: {Sector}")
        if not local:
            base_url = "https://archive.stsci.edu/missions/tess/engineering"
            filelist = [
                file for file in self.listFD(base_url, f"sector{Sector:02d}-{eng}.fits")
            ]
            # TODO Later should potentially code some logic to use the largest number/latest processing
            EngData = fits.open(filelist[0])
        else:
            EngData = self.get_eng_data_local(Sector, eng)
        return EngData

    def Bin_Quat_Cadence(self, cadence, CameraData, Camera):
        """Bin engineering quaternions to observed cadences
        cadence - observed cadences / array of observation midpoints
        CameraData - Quaternion Data For Camera from the quat engineering file
        Camera - Camera Number
        """
        ExpTime = np.double(np.median(np.abs(np.diff(cadence))))
        log.debug(f"Bin_Quat_Cadence: ExpTime={ExpTime}")
        SubArr_Len = int(ExpTime * (60.0 * 60.0 * 24.0))  # convert days to seconds,
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
        # Debug Output
        log.debug(f"Bin_Quat_Cadence: Index Time Length: {Time_Len}")
        log.debug(f"Bin_Quat_Cadence: Source Time Length: {len(CameraData)}")

        for i in range(Time_Len - 1):
            log.debug(
                f"Bin_Quat_Cadence: min: {min_ind} max: {max_ind} max_len: {Time_Len} SubArr_Len: {SubArr_Len}"
            )
            log.debug(f"i: {i} ")

            if max_ind > (Source_Len - 1):
                max_ind = Source_Len - 1
                log.debug(f"\tAdjusting Max Ind: {max_ind}")

            t0 = CameraData["Time"][min_ind]
            t1 = CameraData["Time"][max_ind]
            if t0 > cadence[i]:
                log.debug(
                    f"Index: {i}/{Time_Len} , min, max_ind / source_len: {min_ind}, {max_ind} / {Source_Len} min, max / Time:{t0}, {t1}/ {cadence[i]}"
                )

            if i < (Time_Len - 1) and (max_ind < Source_Len - 1):
                while CameraData["Time"][max_ind] < cadence[i + 1]:
                    max_ind = max_ind + SubArr_Len
                    log.debug(f"\tAdjusting Max Ind: {max_ind}")

                    if max_ind > (Source_Len - 1):
                        max_ind = Source_Len - 1
                        log.debug(f"\tAdjusting Max Ind: {max_ind}")

            if (i > 1) and (min_ind > 0):
                while CameraData["Time"][min_ind] > cadence[i - 1]:
                    min_ind = min_ind - 1
                    log.debug(f"\tAdjusting Min Ind: {min_ind}")

                    if min_ind < 0:
                        min_ind = 0
                        log.debug(f"\tAdjusting Min Ind: {min_ind}")

            if max_ind > (Source_Len - 1):
                max_ind = Source_Len - 1
                log.debug(f"\tAdjusting Max Ind: {max_ind}")

            if min_ind < 0:
                min_ind = 0
                log.debug(f"\tAdjusting Min Ind: {min_ind}")

            SubArr = CameraData[min_ind:max_ind]

            mask = np.nonzero(
                abs(np.double(SubArr["Time"] - cadence[i]))
                < np.double(np.double(0.5) * ExpTime)
            )

            if len(mask[0]) == 0:
                Flag_ZeroMask = Flag_ZeroMask + 1

            # Debug Output
            log.debug(
                f"\t\t\t\tBin_Quat_Cadence: Warning: 0 Quaternions found in the time bin!"
            )
            log.debug(f"Bin_Quat_Cadence: Min Time : {CameraData['Time'][min_ind]}")
            log.debug(f"Bin_Quat_Cadence: Max Time: {CameraData['Time'][max_ind]}")
            log.debug(f"Bin_Quat_Cadence: Cadence: {cadence[i]}")
            log.debug(f"Bin_Quat_Cadence: ExpTime: {ExpTime}")
            log.debug(
                f"Bin_Quat_Cadence: ArrEq: {abs(np.double(SubArr['Time'] - cadence[i])) < np.double(np.double(0.5) * ExpTime)}"
            )
            log.debug(f"Bin_Quat_Cadence: DMin: {cadence[i] - 0.5 * ExpTime}")
            log.debug(f"Bin_Quat_Cadence: DMax: {cadence[i] + 0.5 * ExpTime}")
            log.debug(f"Bin_Quat_Cadence: SubArr: {SubArr['Time']}")
            log.debug(
                f"Bin_Quat_Cadence: Abs(SubArr-Cadence): {abs(np.double(SubArr['Time'] - cadence[i]))}"
            )
            log.debug(f"Bin_Quat_Cadence: Mask: {mask}")
            log.debug(f"Bin_Quat_Cadence: Mask_Len: {len(mask[0])}")
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
                BinnedQuats[7][i] = sigma_clipped_stats(
                    SubArr[f"C{Camera:d}_Q1"][mask]
                )[2]

                BinnedQuats[8][i] = np.median(SubArr[f"C{Camera:d}_Q2"][mask])
                BinnedQuats[9][i] = np.std(SubArr[f"C{Camera:d}_Q2"][mask])
                BinnedQuats[10][i] = sigma_clipped_stats(
                    SubArr[f"C{Camera:d}_Q2"][mask]
                )[2]

                BinnedQuats[11][i] = np.median(SubArr[f"C{Camera:d}_Q3"][mask])
                BinnedQuats[12][i] = np.std(SubArr[f"C{Camera:d}_Q3"][mask])
                BinnedQuats[13][i] = sigma_clipped_stats(
                    SubArr[f"C{Camera:d}_Q3"][mask]
                )[2]

                BinnedQuats[14][i] = np.median(SubArr[f"C{Camera:d}_Q4"][mask])
                BinnedQuats[15][i] = np.std(SubArr[f"C{Camera:d}_Q4"][mask])
                BinnedQuats[16][i] = sigma_clipped_stats(
                    SubArr[f"C{Camera:d}_Q4"][mask]
                )[2]

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
            # print(f"i/Time_Len: {i}/{Time_Len} min_ind, max_ind / source_leb {min_ind}, {max_ind} / {Source_Len}")
            max_ind = int(min_ind) + int(SubArr_Len)
            if max_ind > (Source_Len - 1):
                max_ind = Source_Len - 1
                if t0 > cadence[i]:
                    log.debug(f"\tAdjusting Max Ind: {max_ind}")

        if Flag_ZeroMask > 0:
            log.info(
                f"\t\t\t\tWARNING: Zero Quaternions found for {Flag_ZeroMask} observed cadences"
            )
        return BinnedQuats

    def EMI_Extension_Dict(self, ext_string):
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

    def Bin_EMI_Cadence(self, cadence, EMIData, Camera, type):
        """Bin the EMI information to the observed cadences
        cadence - observed time midpoint array / end-product observation cadences
        EMIData - EMI data file from the emi engineering file
        Camera - camera number
        type - observing cadence for cadence
        """
        ExpTime = np.double(
            abs(cadence[1] - cadence[0])
        )  # Should this be hard-coded instead? prbly med
        if isinstance(type, int):
            type = self.typedict[type]
        match type:
            case "020":
                SubArr_Len = int(4)
            case "120":
                SubArr_Len = int(4)
            case "FFI":
                SubArr_Len = int(
                    (
                        ExpTime
                        / (EMIData[1].data["TIME"][1] - EMIData[1].data["TIME"][0])
                    )
                    * 4
                )
            case _:
                log.info("Error - Bin_EMI_Cadence can't match cadence type")
                # throw a logging warning instead

        Time_Len = len(cadence)
        Source_Len = len(EMIData[1].data["Time"])
        BinnedEMI = np.zeros(
            (10, Time_Len), dtype=np.double
        )  # ( (Earth + Moon) ** (Distance + Elev. + Azim.) + SC(Earth + Moon)*(Ele + Az)
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
                BinnedEMI[0][i] = EMIData[
                    self.EMI_Extension_Dict("emi.earth_distance")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Earth Camera Elev
                BinnedEMI[1][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.cam{Camera}_earth_boresight_angle")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # EarthCamera Azim
                BinnedEMI[2][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.cam{Camera}_earth_azimuth")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Moon Distance
                BinnedEMI[3][i] = EMIData[
                    self.EMI_Extension_Dict("emi.moon_distance")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Moon Camera Elev
                BinnedEMI[4][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.cam{Camera}_moon_boresight_angle")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Moon Camera Azim
                BinnedEMI[5][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.cam{Camera}_moon_azimuth")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Earth Elevation Above Spacecraft
                BinnedEMI[6][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.earth_sc_elevation")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Earth Azimuthal orientation Spacecraft
                BinnedEMI[7][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.earth_sc_azimuth")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Moon Spacecraft Elevation
                BinnedEMI[8][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.moon_sc_elevation")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
                # Moon Azimut Elevation
                BinnedEMI[9][i] = EMIData[
                    self.EMI_Extension_Dict(f"emi.earth_sc_azimuth")
                ].data[min_ind:max_ind]["VALUE"][match_ind]
            else:
                log.info("Bin_EMI_Cadence can't match a datapoint")
                # throw a real logging.warning

            min_ind = min_ind + int(match_ind)
            max_ind = int(min_ind) + int(SubArr_Len)
            if max_ind > Source_Len:
                max_ind = Source_Len - 1

        return BinnedEMI

    def Bin_Quat_Camera(self, CameraData, TimeArr, Camera):
        """Bin quaternions for all cadences for a given camera
        CameraData - quaternion data from the quat engineering file
        TimeArr - observed cadence midpoints from tpfs
        Camera - Camera Number
        """
        Bin20 = None
        if TimeArr[0] is not None:
            log.info("\t\tBinning 20s Quaternion Data")
            Bin20 = self.Bin_Quat_Cadence(TimeArr[0][0], CameraData, Camera)
        log.info("\t\tBinning 120s Quaternion Data")
        Bin120 = self.Bin_Quat_Cadence(TimeArr[1][0], CameraData, Camera)
        log.info("\t\tBinning FFI Quaternion Data")
        BinFFI = self.Bin_Quat_Cadence(TimeArr[2][0], CameraData, Camera)

        return Bin20, Bin120, BinFFI

    def Bin_EMI_Camera(self, EMIData, TimeArr, Camera):
        """Bin emi info for all cadences for a given camera
        EMIData - EMI data from the emi engineering file
        TimeArr - observed cadence midpoints from tpfs
        Camera - Camera Number
        """
        EMI20 = None
        if TimeArr[0] is not None:
            log.info("\t\tBinning 20s EMI Data")
            EMI20 = self.Bin_EMI_Cadence(TimeArr[0][0], EMIData, Camera, "020")

        log.info("\t\tBinning 120s EMI Data")
        EMI120 = self.Bin_EMI_Cadence(TimeArr[1][0], EMIData, Camera, "120")

        log.info("\t\tBinning FFI EMI Data")
        EMIFFI = self.Bin_EMI_Cadence(TimeArr[2][0], EMIData, Camera, "FFI")

        return EMI20, EMI120, EMIFFI

    def get_break_times(self, times):
        break_inds = np.where(np.diff(np.asarray(times)) > 1 / 24)[0] + 1
        break_times = np.asarray(
            [np.mean(times[ind - 1 : ind + 1]) for ind in break_inds]
        )
        return break_times

    def get_segment_label(self, times, break_times):
        """Create labels for segments - using the quality flag and timing array

        right now, an 'segment' is every stretch of observations between earth point
        quality flags (=8), with each earth point stretch causing an increment in
        segment label
        """
        label = 1
        # cut should not be nescessary
        cut_times = times[times != 0]
        segment_labels = np.zeros_like(cut_times)

        for btime in break_times:
            if label == 1:
                segment_labels[cut_times < btime] = label
            label += 1
            segment_labels[cut_times > btime] = label
        if np.any(segment_labels == 0):
            raise RuntimeWarning(
                "get_segment_label warning: Not all times have a segment label"
            )

        return segment_labels

    def write_file_header(self, fname, Sector, Camera, Cadence):
        with open(fname, "w") as f:
            f.write(f"# TESS Quaternions downsampled to end-user cadences\n")
            f.write(f"# Sector: \t{Sector}\n")
            f.write(f"# Camera: \t{Camera}\n")
            f.write(f"# Cadence:\t{Cadence}\n")
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
                f"# MidTime: The CCD1 exposure midpoint in spacecraft time (e.g. 'TIME' - 'TIMECORR from a SPOC TPF)'. Other CCDs will have small (~0.2-1)s offsets from this due this due to sequential reads.\n"
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
            f.write(f"# Earth_Distance: Distance to Earth in Earth Radii \n")
            f.write(
                f"# Earth_Camera_Angle: Angle of Earth from Camera Boresight in Degrees \n"
            )
            f.write(
                f"# Earth_Camera_Azimuth: Azimuth of Earth around Camera Boresight in Degrees \n"
            )
            f.write(f"# Moon_Distance: Distance to Moon in Earth Radii \n")
            f.write(
                f"# Moon_Camera_Angle: Angle of Moon from Camera Boresight in Degrees \n"
            )
            f.write(
                f"# Moon_Camera_Azimuth: Azimuth of Moon around Camera Boresight in Degrees \n"
            )
            f.write(
                f"# Earth_Spacecraft_Angle: Angle of Earth from Spacecraft Boresight in Degrees \n"
            )
            f.write(
                f"# Earth_Spacecraft_Azimuth: Azimuth of Earth around Spacecraft Boresight in Degrees \n"
            )
            f.write(
                f"# Moon_Spacecraft_Angle: Angle of Moon from Sacecraft Boresight in Degrees \n"
            )
            f.write(
                f"# Moon_Spacecraft_Azimuth: Azimuth of Moon around Spacecraft Boresight in Degrees \n\n"
            )

            if Cadence == self.typedict[3]:
                f.write(
                    f"# FFI_File: The FFI file assosciated with this observing cadence\n"
                )
            f.write(f"# Processing Date-Time: {date.today()}\n\n")

    def write_vector_sector_camera(
        self, BinnedQuats, BinnedEMI, FFIList, TimeArr, Sector, Camera
    ):
        """Write TESSVectors Info to CSV files for all observed cadences"""
        # Binned_Dir = "TESSVectors_products"
        for Quat, EMI, Time, i in zip(BinnedQuats, BinnedEMI, TimeArr, [1, 2, 3]):
            fname_root = self._vector_file(self.typedict[i], Sector, Camera).split(".")[
                0
            ]
            fname = fname_root + ".xz"
            log.info(
                f"\t\tWriting Sector: {Sector} Camera: {Camera} Cadence {self.typedict[i]} to:"
            )
            log.info(f"\t\t\t to {fname}")

            if Time is not None:
                # First write the File Header
                # this will over-write by default

                self.write_file_header(fname, Sector, Camera, self.typedict[i])
                df = pd.DataFrame(
                    data={
                        "Cadence": Time[2],
                        "MidTime": Time[0],
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
                        "Earth_Spacecraft_Angle": EMI[6],
                        "Earth_Spacecraft_Azimuth": EMI[7],
                        "Moon_Spacecraft_Angle": EMI[8],
                        "Moon_Spacecraft_Azimuth": EMI[9],
                    }
                )
                if i == 3:
                    df["FFIFile"] = FFIList

                df.to_csv(fname, index=False, mode="a", compression=None)

                df = pd.read_csv(fname, comment="#", index_col=False, compression=None)
                df.to_csv(fname, index=False, mode="w", compression="xz")

    def create_vectors_sector(self, Sector, check_exists=True):
        """For a given sector, create TESSVectors Information and write to a CSV file"""
        # assume that files don't exist, unless we're checking.  could rename this overwrite?
        files_exist = False
        cutoff_20s = 27
        if check_exists & self.check_exists:
            Camera = [1, 2, 3, 4]

            if Sector >= cutoff_20s:
                i = [1, 2, 3]
            else:
                i = [2, 3]

            files_exist = True
            for item in product(i, np.atleast_1d(Camera)):
                file_check = os.path.isfile(
                    self._vector_file(self.typedict[item[0]], Sector, item[1])+'.xz'
                )
                files_exist = files_exist and file_check

        if not files_exist:
            log.info(f"Starting Sector: {Sector}")
            QuatData = self.get_eng_data(Sector, "quat")
            EMIData = self.get_eng_data(Sector, "emi")

            for Camera in [1, 2, 3, 4]:
                camera_files_exist = True
                if check_exists & self.check_exists:
                    camera_files_exist = True
                    for item in product(i, np.atleast_1d(Camera)):
                        camera_file_check = os.path.isfile(
                            self._vector_file(self.typedict[item[0]], Sector, item[1])+'.xz'
                        )

                        camera_files_exist = camera_files_exist and camera_file_check
                if not camera_files_exist:
                    log.info(f"\tStarting Camera: {Camera}")
                    TimeArr, FFIList = self.get_camera_sector_cadences(Sector, Camera)
                    BinnedQuats = self.Bin_Quat_Camera(
                        QuatData[Camera].data, TimeArr, Camera
                    )
                    BinnedEMI = self.Bin_EMI_Camera(EMIData, TimeArr, Camera)
                    self.write_vector_sector_camera(
                        BinnedQuats, BinnedEMI, FFIList, TimeArr, Sector, Camera
                    )
                else:
                    log.info(f"Sector: {Sector} Camera: {Camera} files exist, skipping")

        else:
            log.info(f"Sector: {Sector} files exist, skipping")
