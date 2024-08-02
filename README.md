# tessvectors 

# PRE-RELEASE 
0.1.0 - Open Beta 

## A repository that takes Transiting Exoplanet Survey Satelitte (TESS) engineering data and transforms it into convenient value-added products for use in lightcurve detrending, image characterising, and more.   

## tessvectors Introduction

The tessvectors repository is an effort to take TESS mission engineering products and translate them into a more convenient, value-added form factor that is convenient for use by end-users.  These files are essentially csv files that can be read by pandas (though are currently being compressed using lzma compression)

The full available vectors and diagnostic files can be found at:

https://heasarc.gsfc.nasa.gov/docs/tess/data/TESSVectors/

We have a [tutorial notebook](docs/tessvectors_tutorial.ipynb) on how to interact with the vectors.  The `tessvectors` package also contains a few convenience functions in a lightweight API that may help you interact with these data, and are demonstrated here. 

We also have a  [demo notebook](docs/tessvectors_demo.ipynb) that demonstrates a few things you can do with the information contained in these files.

A demonstration of products can be found in the [Products](Products/) folder.  

## tessvectors usage
tessvectors vector files can be downloaded from https://heasarc.gsfc.nasa.gov/docs/tess/data/TESSVectors/ and read in using `pandas.read_csv`.  The `tessvectors` package can be installed using:

```
pip install git+https://github.com/tessgi/tessvectors.git
```

A vector file can be read using the `tessvectors.get_vectors` function and supplying a fits file, path, or specific cadence/sector/camera information:

```
import tessvectors
vector = tessvectors.getvector(("FFI", 1, 4))
```

We also have a few convenience functions.  One involves using a given time, Sector, and Camera to look up a specific FFI:

```
ffi_loc = tessvectors.getffi((1337.370729227238,1,4,1))
```

Please see the [tutorial notebook](docs/tessvectors_tutorial.ipynb) for more details.  

## tessvectors Information
The information contained inside of these files is:

    - Cadence #: Cadence index from the source tpf
    - MidTime: The exposure midpoint in spacecraft time without barycentric correction (i.e. tpf.time - tpf.timecorr)
        - This has been benchmarked to CCD1.  Due to staggered read, other CCDs midpoints might differ by up to ~1s
    - ExpTime: The final cadence binning (20s/120s/FFI)
    - Sector: The TESS observing Sector for the source data
    - Camera: The TESS camera for the source data
    - Quat_Start: The timestamp of the earliest quaternion used in the bin
    - Quat_Stop: The timestamp of the last quaternion used in the bin
    - Quat_MIN_FOM: The worst Figure of Merit from the source quaternions
    - Quat_MIN_NUM_GSUSED: The lowest number of guide stars used in the source quaternions
    - Quat_NBinned: The number of quaternions binned into this final result.
    - Quat[1-4]_Med: The Quaternion #[1-4] median value from the binned values 
    - Quat[1-4]_StdDev: The standard deviation of Quaternion #[1-4] binned values
    - Quat[1-4]_SigClip: The Sigma-Clipped Standard Deviation of Quaternion #[1-4] binned values
    - Earth_Distance: Distance to Earth in Earth Radii
    - Earth_Camera_Angle: Angle of Earth from Camera Boresight in Degrees
    - Earth_Camera_Azimuth: Azimuth of Earth around Camera Boresight in Degrees
    - Moon_Distance: Distance to Moon in Earth Radii
    - Moon_Camera_Angle: Angle of Moon from Camera Boresight in Degrees
    - Moon_Camera_Azimuth: Azimuth of Moon around Camera Boresight in Degrees
    - (FFI Cadences Only) FFIFile: The FFI file assosciated with this observing cadence 

The information contained in these files generally comes from quaternion (\*-quat.fits) and earth-moon information (\*-emi.fits) [engineering files](https://archive.stsci.edu/missions-and-data/tess/data-products.html#mod_eng) that have been reprocessed to present results at the same time-cadence/binning as end usser TPFs'/lightcurves.

There is one TESS-Vectors file for every Cadence(20-second/120-second/FFI) and Camera/Sector.  We have also created some diagnostic plots to help guide the understanding of TESS observations for a given sector.  

### tessvectors Data Usage Notes
While we have endeavoured to make these files as simple and straightforward to use as possible, there are a few points of note that one should keep in mind when using these files: 

    - MidTime: tessvectors data products are created on a per-camera basis as that is the referance frame for the engineering files, however the CCD readouts can have small variations (<1s between CCD's 1-4) and so may not line up exactly with these values.  The cadence numbers should align with cadence numbers from mission TPFs and can be used to cross-index data as well. Currently CCD 1 for all cameras is being used as a benchmark.
    - Not all camera's will have quaternions at all times.  NaN's have been introduced where this data is unavailable.  In this event other cameras are used for guiding (typically cameras 1 & 4 will be used to guide TESS).

### tessvectors processing
The functionality to create the tessvectors files may be found in the [processing/](src/tessvectors/processing/) subpackage.  Please see the processing [README](src/tessvectors/processing/README.md) for further information.  

### Credits

Initial development done by [Tyler Pritchard](https://github.com/tylerapritchard), [Christina Hedges](https://github.com/christinahedges) the [TESS Science Support Center](https://heasarc.gsfc.nasa.gov/docs/tess/), and the [MIT TESS Science Operations Team](https://tess.mit.edu/). 
