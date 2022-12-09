# EarthquakeFrequency

This code gets earthquake information using python wrapper ‘libcomcat’ 
developed by the USGS (available at https://github.com/usgs/libcomcat) 
and download ShakeMap information.

The number of stations recorded the earthquake can be filtered based on
distance from epicenter and shaking intensity at site. The seismograms 
from these stations are downloaded via ObsPy developed by IRIS 
(available at https://github.com/obspy/obspy). Both instrument response
and baseline correction is automatically done to the downloaded seismograms.

To get instantenous frequency at the PGA of these seismograms, the code applies
Hilbert-Huang transform method (HHT) method (https://en.wikipedia.org/wiki/Hilbert–Huang_transform)
to the 10 secs before and after the PGA.

The HHT works by first decomposing a signal into a series of intrinsic mode functions (IMFs) 
using the EMD method. This process is called sifting which was first developed by 
Huang et al. (1998). We use python module ‘pyEMD’ 
(available at https://pyemd.readthedocs.io/en/latest/index.html) for computing this process 
which utilizes the algorithm developed by Rilling et al. (2003).
