import tessvectors
from tessvectors.processing import makevectors
import pandas as pd
import lightkurve as lk
import os


def test_getvector():
    # can we get a vector via specifying a Cadence, Sector Camera
    df = tessvectors.getvector(("FFI", 1, 4))
    assert isinstance(df, pd.DataFrame)

    tpf = lk.search_targetpixelfile(
        "Pi Men C", sector=1, exptime="short", author="SPOC"
    ).download()

    # cand we get a vevtor vid passing a hdu?
    df = tessvectors.getvector(tpf.hdu)
    assert isinstance(df, pd.DataFrame)

    # can we get a vector via specifying a file?
    df = tessvectors.getvector(tpf.path)
    assert isinstance(df, pd.DataFrame)

    # can we download, and read-in a vector file
    df = tessvectors.getvector(("FFI", 1, 4), download=True)
    assert os.path.isfile(makevectors()._vector_base_file("FFI", 1, 4) + ".csv")

def test_getffi():
    ffi_loc = tessvectors.getffi((1337.370729227238,1,4,1))
    assert isinstance(ffi_loc, str)
