README file for TESSVectors:Value-added engineering products for use detrending, image characterisation, and more.    
MAST webpage: 
Refer to this HLSP with DOI:

# TESSVectors Introduction

The TESS-Vectors repository is an effort to take TESS mission engineering products and translate them into a more convenient, value-added form factor that is convenient for use by end-users.  

The information contained in these files primarily comes from quaternion (\*-quat.fits) and earth-moon information (\*-emi.fits) [engineering files](https://archive.stsci.edu/missions-and-data/tess/data-products.html#mod_eng) that have been reprocessed to present results at the same time-cadence/binning as end usser TPFs'/lightcurves.



# Data Products

There is one TESS-Vectors file for every Cadence(20-second/120-second/FFI) and Camera/Sector.  We have also created some diagnostic plots to help guide the understanding of TESS observations for a given sector.  

The naming schema for the TESSVectos csv files are: 
    
    - TessVectors_S{Sector}_C{Camera}_{Observing Cadence}.csv
        - Sector = zero-padded four-digit sector number, e.g. 0001
        - Camera = Integer Camera Number (1-4)
        - Observing Cadence = three letter string 020 / 120 / FFI corresponding to the 20-second, 120-second, and FFI TESS Observing Cadences

The information contained inside of these files is:

    - Cadence #: Cadence index from the source tpf
    - MidTime: The exposure midpoint in spacecraft time (i.e. tpf.time - tpf.timecorr)
    - TimeCorr: The Time Correction for spacecraft to Barycentric time at that cadence
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


# Source Code
The source code that was used to create these files can be found at: 

# Contributors
Initial development done by [Tyler Pritchard](https://github.com/tylerapritchard), [Christina Hedges](https://github.com/christinahedges) the [TESS Science Support Center](https://heasarc.gsfc.nasa.gov/docs/tess/), and the [MIT TESS Science Operations Team](https://tess.mit.edu/)