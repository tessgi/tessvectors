import tessvectors
import pandas as pd
import lightkurve as lk
import os


def test_vector():
    # can we get a vector via specifying a Cadence, Sector Camera
    df = tessvectors.vector(("FFI", 1, 4))
    assert isinstance(df, pd.DataFrame)

    tpf = lk.search_targetpixelfile(
        "Pi Men C", sector=1, exptime="short", author="SPOC"
    ).download()

    # cand we get a vevtor vid passing a hdu?
    df = tessvectors.vector(tpf.hdu)
    assert isinstance(df, pd.DataFrame)

    # can we get a vector via specifying a file?
    df = tessvectors.vector(tpf.path)
    assert isinstance(df, pd.DataFrame)

    # can we download, and read-in a vector file
    df = tessvectors.vector(("FFI", 1, 4), download=True)
    assert os.path.isfile(
        tessvectors.processing()._vector_base_file("FFI", 1, 4) + ".csv"
    )
