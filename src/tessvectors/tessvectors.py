import astropy.io.fits as fits
import pandas as pd
import logging
import sys

from multiprocessing.pool import Pool
import multiprocessing as mp

from .processing import processing
from .diagnostics import diagnostics

logging.basicConfig(level=logging.INFO, stream=sys.stdout)
log = logging.getLogger("tessvectors")
log.setLevel(logging.INFO)


def _get_Cadence_Sector_Camera(hdu):
    Sector = hdu[0].header["SECTOR"]
    Camera = hdu[0].header["CAMERA"]
    t = hdu[1].data["TIME"]
    dt = (t[1] - t[0]) * 24 * 60 * 60
    if (dt < 30) and (dt > 19):
        Cadence = "020"
    if (dt > 119) and (dt < 135):
        Cadence = "120"
    if dt > 199:
        Cadence = "FFI"
    return Cadence, Sector, Camera


def _resolve_remote_file(Cadence, Sector, Camera):
    remote_base = (
        "https://github.com/tylerapritchard/TESSVectors/raw/main/Products/Vectors/"
    )
    cadence_folder = f"{Cadence}_Cadence/"
    file_base = f"TessVectors_S{Sector:03d}_C{Camera}_{Cadence}.xz"
    fname = remote_base + cadence_folder + file_base
    return fname


def vector(arg, Cadence=None, Sector=None, Camera=None, download=False):
    if isinstance(arg, tuple):
        # Assume a tupe is a set of Cadence/Sector/Camera
        if len(arg) == 3:
            Cadence, Sector, Camera = arg
        else:
            raise KeyError(
                "Must Provide (Cadence(str), Sector(int), Camera(int) when supplying a tuple)"
            )
    elif isinstance(arg, str):
        # Assume a string is a file and read it
        with fits.open(arg) as hdu:
            Cadence, Sector, Camera = _get_Cadence_Sector_Camera(hdu)
    # if we are given a hdulist object, assume this is a lc/tpf and determine the
    # Cadence/Sector/Camera
    elif isinstance(arg, fits.hdu.HDUList):
        Cadence, Sector, Camera = _get_Cadence_Sector_Camera(arg)
    else:
        raise TypeError("Unknown Input Supplied")

    fname = _resolve_remote_file(Cadence, Sector, Camera)
    logging.info(f"Retrieving TESSVector from: {fname}")
    vector = pd.read_csv(fname, comment="#", index_col=False)

    if download:
        output_file = processing()._vector_base_file(Cadence, Sector, Camera)
        # lets assume people want uncompressed csv files by default
        output_file = f"{output_file}.csv"
        processing().write_file_header(output_file, Sector, Camera, Cadence)
        vector.to_csv(output_file, mode="a", compression=None)

    return vector


### Multiprocessing & convenience functions
def process_sector(Sector):
    try:
        processing().create_vectors_sector(Sector)
    except:
        log.warning(
            f"\t\t\t Warning, Creating Vector for Sector: {Sector} Failed!"
        )  # add a real warning
        pass
    try:
        diagnostics().create_diagnostics_sector(Sector)
    except:
        log.warning("\t\t\t Warning, Plotting Failed!")  # add a real warning
        pass
    return (Sector, True)


def run_bulk_vectors(sector_min=1, sector_max=78, processes=4):
    if not processes:
        processes = 7
    sector_list = range(sector_min, sector_max)
    pool = Pool(processes=processes)
    res = []
    for result in pool.imap_unordered(
        processing().create_vectors_sector, sector_list, chunksize=1
    ):
        res = [res, result]
    pool.close()


def run_bulk_processing(sector_min=1, sector_max=77):
    sector_list = range(sector_min, sector_max)
    pool = Pool(processes=mp.cpu_count())
    res = []
    for result in pool.imap_unordered(process_sector, sector_list, chunksize=1):
        res = [res, result]
    pool.close()


def run_bulk_diagnostics(self, sector_min=1, sector_max=78, processes=7):
    if not processes:
        processes = 7
    sector_list = range(sector_min, sector_max)
    pool = Pool(processes=processes)
    res = []
    for result in pool.map(diagnostics().create_diagnostics_sector, sector_list):
        res = [res, result]
    pool.close()