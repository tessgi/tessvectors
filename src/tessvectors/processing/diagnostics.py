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

import numpy as np

import pandas as pd

import lightkurve as lk
from lightkurve import TessQualityFlags

from . import makevectors

mpl.rcParams["agg.path.chunksize"] = 10000000

logging.basicConfig(level=logging.INFO, stream=sys.stdout)
log = logging.getLogger("tessvectors")
log.setLevel(logging.INFO)


class diagnostics(object):
    def __init__(
        self,
        TESSVectors_Products_Base="Products",
        TESSVectos_local=False,
        TESSVectors_Local_ENG_Path="Eng",
        TESSVectors_Local_tpf_fast_path="SourceData/fast",
        TESSVectors_Local_tpf_short_path="SourceData/short/",
        TESSVectors_Local_tpf_ffi_path="SourceData/FFI/",
        check_exists=True,
        compression=False,
    ):
        self.TESSVectors_Products_Base = TESSVectors_Products_Base
        self.TESSVectos_local = TESSVectos_local
        self.TESSVectors_Local_ENG_Path = TESSVectors_Local_ENG_Path
        self.TESSVectors_Local_tpf_fast_path = TESSVectors_Local_tpf_fast_path
        self.TESSVectors_Local_tpf_short_path = TESSVectors_Local_tpf_short_path
        self.TESSVectors_Local_tpf_ffi_path = TESSVectors_Local_tpf_ffi_path
        self.check_exists = check_exists
        self.typedict = {1: "020", 2: "120", 3: "FFI"}
        self.compression = compression

        if not compression:
            self.vector_extenstion = ".csv"
        else:
            self.vector_extenstion = ".xz"

        makevectors().make_dir_structure()

    def create_diagnostics_sector(self, Sector):
        # CameraCadence
        # Sector, Camera, Cadence = SectorCameraCadence
        log.info(f"Creating Diagnostics for Sector: {Sector}")
        for Camera in [1, 2, 3, 4]:
            try:
                self.create_diagnostic_periodogram(Sector, Camera)
                for Cadence in [1, 2, 3]:
                    self.create_diagnostic_timeseries(Sector, Camera, Cadence)
                    # Should I create the periodograms from the "raw" 2s data?  probably?
                    self.create_diagnostic_emi(Sector, Camera, Cadence)

            except:
                print(f"Sector {Sector} Diagnostics Failed")

    def _quality_to_color(qual):
        default = TessQualityFlags.create_quality_mask(qual, bitmask="default")
        hard = TessQualityFlags.create_quality_mask(qual, bitmask="hard")
        # hardest = TessQualityFlags.create_quality_mask(qual, bitmask='hardest')

        carr = np.zeros(len(qual))
        # carr[~ hardest] = 0.33
        carr[~hard] = 0.5
        carr[~default] = 1

        return carr

    def _plot_quat(self, axs, time, quat, dev, qual, QuatLabel):
        nsigma = 1

        # norm = np.nanmedian(np.double(quat[~np.isfinite(quat)]))
        norm = np.nanmedian(quat)
        norm_quat = quat / norm
        norm_dev = dev / norm
        # mean_norm_dev = np.nanmedian(dev[~np.isfinite(dev)] / norm)
        mean_norm_dev = np.nanmedian(dev / norm)
        tmin = np.min(time)
        tmax = np.max(time)

        ymin = np.nanmedian(norm_quat) - 3 * mean_norm_dev
        ymax = np.nanmedian(norm_quat) + 3 * mean_norm_dev
        if qual is not None:
            im = axs.imshow(
                np.vstack((qual,)),
                extent=(tmin, tmax, ymin, ymax),
                interpolation="nearest",
                aspect="auto",
                cmap=cm.PuRd,
                vmax=1,
            )
        else:
            im = None

        axs.scatter(
            time,
            norm_quat,
            s=3,
            label="Median Quaternion Value",
            color="k",
            rasterized=True,
        )

        axs.fill_between(
            time,
            norm_quat - (nsigma * norm_dev),
            norm_quat + (nsigma * norm_dev),
            alpha=0.6,
            color="grey",
            rasterized=True,
        )

        axs.set_xlim(tmin, tmax)
        axs.set_ylim(ymin, ymax)

        axs.tick_params(axis="x", labelsize=1)

        axs.set_ylabel(f"{QuatLabel}", weight="bold", size=24)
        axs.set_yticks([])
        axs.set_rasterized(True)
        return im

    def _diagnostic_timeseries_file(self, cadence_name, Sector, Camera):
        return f"{self.TESSVectors_Products_Base}/Diagnostics/Vectors/TessVectors_S{Sector:03d}_C{Camera}_{cadence_name}_Quat.pdf"

    def _diagnostic_power_file(self, Sector, Camera):
        return f"{self.TESSVectors_Products_Base}/Diagnostics/Periodograms/TessVectors_S{Sector:03d}_C{Camera}_QuatPower.pdf"

    def create_diagnostic_timeseries(self, Sector, Camera, Cadence):
        # TODO add segment labels to plot?

        if type(Cadence) != str:
            cadence_name = self.typedict[Cadence]
        else:
            cadence_name = Cadence

        fname = (
            makevectors()._vector_file(cadence_name, Sector, Camera)
            + self.vector_extenstion
        )

        nplots = 3
        if os.path.isfile(fname):
            quatdf = pd.read_csv(fname, comment="#", index_col=False)
            # qual_im = self._quality_to_color(quatdf.Quality)
            qual_im = None

            fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))
            im = self._plot_quat(
                axs[0],
                quatdf.MidTime,
                quatdf.Quat1_Med,
                quatdf.Quat1_StdDev,
                qual_im,
                "Quaternion 1",
            )

            im = self._plot_quat(
                axs[1],
                quatdf.MidTime,
                quatdf.Quat2_Med,
                quatdf.Quat2_StdDev,
                qual_im,
                "Quaternion 2",
            )

            im = self._plot_quat(
                axs[2],
                quatdf.MidTime,
                quatdf.Quat3_Med,
                quatdf.Quat3_StdDev,
                qual_im,
                "Quaternion 3",
            )

            plt.subplots_adjust(hspace=0)
            axs[0].set_title(
                f"TESS Sector {Sector} Camera {Camera} Quaternions",
                weight="bold",
                size=26,
            )
            axs[-1].tick_params(axis="x", labelsize=18)
            axs[-1].set_xlabel("TESS BTJD", weight="bold", size=24)

            # cax = plt.axes([0.92, 0.11, 0.075, 0.77])
            # cbar = plt.colorbar(mappable=im, cax=cax, ticks=[0, 0.5, 1])
            # cbar.ax.set_yticklabels(["Unflagged", "Aggressive", "Conservative"], size=18)
            # cbar.set_label("Data Flagging Level (Lower Is Better)", size=24, weight="bold")

            # plt.tight_layout()
            fout = self._diagnostic_timeseries_file(cadence_name, Sector, Camera)
            # plt.show()
            plt.savefig(fout, dpi=100, bbox_inches="tight")  # , rasterize=True)
            plt.show()
            plt.close(fig)
        else:
            log.warning(f"WARNING: Vector file {fname} does NOT exist")

    def plot_lsperiodogram(self, ax, time, median, std, QuatLabel):
        lc = lk.LightCurve(data={"time": time, "flux": std})
        ls = lc.to_periodogram(maximum_period=13.5)
        ls.plot(ax=ax, lw=0.1, color="k", ylabel=" ", rasterized=True)

        ax.set_ylabel(f"{QuatLabel} Power", weight="bold", size=24)
        # ax.set_yticks([])
        ax.tick_params(axis="y", labelsize=18)
        ax.set_yscale("log")
        ax.set_xscale("log")

        ax.set_ylim(min(ls.power), max(ls.power))
        return

    def create_diagnostic_periodogram(self, Sector, Camera):
        cutoff_20s = 27
        if Sector < cutoff_20s:
            Cadence = 2
        else:
            Cadence = 1

        if type(Cadence) != str:
            cadence_name = self.typedict[Cadence]
        else:
            cadence_name = Cadence
        fname = (
            makevectors()._vector_file(cadence_name, Sector, Camera)
            + self.vector_extenstion
        )
        nplots = 3
        if os.path.isfile(fname):
            if not self.compression:
                quatdf = pd.read_csv(
                    fname, comment="#", index_col=False, low_memory=False
                )
            else:
                quatdf = pd.read_csv(
                    fname,
                    comment="#",
                    index_col=False,
                    compression="xz",
                    low_memory=False,
                )
            fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))
            self.plot_lsperiodogram(
                axs[0],
                quatdf.MidTime,
                quatdf.Quat1_Med,
                quatdf.Quat1_StdDev,
                "Quaternion 1",
            )
            self.plot_lsperiodogram(
                axs[1],
                quatdf.MidTime,
                quatdf.Quat2_Med,
                quatdf.Quat2_StdDev,
                "Quaternion 2",
            )
            self.plot_lsperiodogram(
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

            fout = self._diagnostic_power_file(Sector, Camera)
            plt.savefig(fout, dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            log.warning(f"WARNING: Vector file {fname} does NOT exist")

    def _diagnostic_emi_file(self, Sector, Camera):
        return f"{self.TESSVectors_Products_Base}/Diagnostics/EMI/TessVectors_S{Sector:03d}_C{Camera}_emi.pdf"

    def create_diagnostic_emi(self, Sector, Camera, Cadence):
        if type(Cadence) != str:
            cadence_name = self.typedict[Cadence]
        else:
            cadence_name = Cadence

        fname = (
            makevectors()._vector_file(cadence_name, Sector, Camera)
            + self.vector_extenstion
        )

        nplots = 3
        if os.path.isfile(fname) and cadence_name == "FFI":
            if not self.compression:
                quatdf = pd.read_csv(
                    fname, comment="#", index_col=False, low_memory=False
                )
            else:
                quatdf = pd.read_csv(
                    fname,
                    comment="#",
                    index_col=False,
                    compression="xz",
                    low_memory=False,
                )

            # qual_im = self._quality_to_color(quatdf.Quality)
            qual_im = None

            fig, axs = plt.subplots(nplots, 1, figsize=(15, nplots * 10))

            axs[0].scatter(
                quatdf.MidTime, quatdf.Earth_Distance, label="Earth", color="seagreen"
            )
            axs[0].scatter(
                quatdf.MidTime,
                quatdf.Moon_Distance,
                label="Moon",
                color="darkturquoise",
            )
            axs[0].legend(
                prop={"weight": "bold", "size": 24}, scatterpoints=3, markerscale=2
            )

            tmin = np.min(quatdf.MidTime)
            tmax = np.max(quatdf.MidTime)

            ymin = min(min(quatdf.Earth_Distance), min(quatdf.Moon_Distance))
            ymax = max(max(quatdf.Earth_Distance), max(quatdf.Moon_Distance))
            if qual_im is not None:
                im = axs[0].imshow(
                    np.vstack((qual_im,)),
                    extent=(tmin, tmax, ymin, ymax),
                    interpolation="nearest",
                    aspect="auto",
                    cmap=cm.PuRd,
                    vmax=1,
                )

            axs[1].scatter(
                quatdf.MidTime,
                quatdf.Earth_Camera_Angle,
                label="Earth",
                color="seagreen",
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

            if qual_im is not None:
                im = axs[1].imshow(
                    np.vstack((qual_im,)),
                    extent=(tmin, tmax, ymin, ymax),
                    interpolation="nearest",
                    aspect="auto",
                    cmap=cm.PuRd,
                    vmax=1,
                )

            axs[2].scatter(
                quatdf.MidTime,
                quatdf.Earth_Camera_Azimuth,
                label="Earth",
                color="seagreen",
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

            ymin = min(
                min(quatdf.Earth_Camera_Azimuth), min(quatdf.Moon_Camera_Azimuth)
            )
            ymax = max(
                max(quatdf.Earth_Camera_Azimuth), max(quatdf.Moon_Camera_Azimuth)
            )

            if qual_im is not None:
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

            # cax = plt.axes([0.92, 0.11, 0.075, 0.77])
            # cbar = plt.colorbar(mappable=im, cax=cax, ticks=[0, 0.5, 1])
            # cbar.ax.set_yticklabels(["Unflagged", "Aggressive", "Conservative"], size=18)
            # cbar.set_label("Data Flagging Level (Lower Is Better)", size=24, weight="bold")

            fout = self._diagnostic_emi_file(Sector, Camera)
            plt.savefig(fout, dpi=300, bbox_inches="tight")
            # plt.show()
            plt.close(fig)
        else:
            log.warning(f"WARNING: Vector file {fname} does NOT exist")
