![alt tag](http://dbuscombe-usgs.github.io/figs/pysesa_colour.jpg)

PySESA - a Python framework for Spatially Explicit Spectral Analysis

PySESA is an open-source project dedicated to provide a generic Python framework
for spatially explicit statistical analyses of point clouds and other geospatial data,
in the spatial and frequency domains, for use in the geosciences

The program is detailed in:
   Buscombe, D. (2016) "spatially explicit spectral analysis of point clouds and geospatial data", computers and geosciences 86, 92-108, 10.1016/j.cageo.2015.10.004.

For more information visit http://dbuscombe-usgs.github.io/pysesa/

    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information,
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only
    and does not imply endorsement by the U.S. government.

### Contributing & Credits

 Author |    Daniel Buscombe
 ------ | ---------------
        |  Grand Canyon Monitoring and Research Center
        | United States Geological Survey
        | Flagstaff, AZ 86001
        | Now at Northern ARizona University. daniel.buscombe@nau.edu
Version: 2   |  Revision: Nov, 2019

For latest code version please visit:

```
https://github.com/dbuscombe-usgs
```

This function is part of pysesa software
This software is in the public domain because it contains materials that originally came
from the United States Geological Survey, an agency of the United States Department of Interior.
For more information, see the official USGS copyright policy at

```
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
```

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government.

### License
GNU Lesser General Public License, Version 3
```
http://www.gnu.org/copyleft/lesser.html
```


## Setup

### Install
Clone the repo from github:

```
git clone git@github.com:dbuscombe-usgs/pysesa.git
```

change directory:
```
cd pysesa
```

Create a conda environment that contains all the libraries necessary to run pysesa
```
conda env create -f conda_env/owg.yml
```

Activate the environment:
```
conda activate pysesa
```

Compile the cython modules in place. Modules are already provided for python 3.7 on 64-bit Windows 10. Otherwise, you will need to recompile

```
python setup.py build_ext --inplace
```

### Test

A test can be carried out by running the supplied script:

```
python -c "import pysesa; pysesa.test()"
```

which carries out the following operations:

```
   # general settings   
   infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_100000pts.xyz'

   out = 0.5 #m output grid
   detrend = 4 #ODR plane
   proctype = 1 #Processing spectral parameters (no smoothing)
   mxpts = 1024 # max pts per window
   res = 0.05 #cm internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # Hann taper
   prc_overlap = 100 # 100% overlap between successive windows
   minpts = 64 # min pts per window

   pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)

```

### Support

Please use github issues

### Contents

The programs in this package are as follows:

1. read: read a 3-column space, comma or tab delimited text file
2. partition: partition a Nx3 point cloud into M windows of nx3 points with specified spacing between centroids of adjacent windows and with specified overlap between windows.
3. detrend: returns detrended amplitudes of a Nx3 point cloud
4. sgolay: returns the Savitsky-Golay digital filter of a 2D signal
5. spatial: calculate spatial statistics of a Nx3 point cloud
6. RunningStats: called by \spatial} to compute sigma, skewness and kurtosis
7. lengthscale: calculates the integral lengthscale of a Nx3 point cloud
8. spectral: calculate spectral statistics of a Nx3 point cloud
9. process: allows control of inputs to all modules (full workflow)
10. write: write program outputs to a comma delimited text file
11. test: program testing suite

These are all command-line/modular programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options


### Read
    '''
    Custom fast (up to 3.5x faster than numpy's genfromtxt) txt file to numpy array
    accepts comma, tab or space delimited files of 3 columns: x, y, and amplitude

    Syntax
    ----------
    pts = pysesa_read.txtread(infile)

    Parameters
    ----------
    infile : str
   	3-column ASCII file containing Nx3 point cloud

    Returns
    ----------
    data: ndarray
   	Nx3 point cloud, 32 bit precision

    '''

### Partition
    '''
    Partition a Nx3 point cloud into M windows of nx3 points
    with specified spacing between centroids of adjacent windows
    and with specified overlap between windows.
    Implemented using a binary search tree for fast nearest neighbour
    point check with boundary pruning

    Syntax
    ----------
    nr_pts = pysesa.partition(toproc, out, res, mxpts, minpts, prc_overlap).getdata()

    Parameters
    ----------
    toproc : ndarray
   	Nx3 point cloud

    Other Parameters
    ----------
    out : float, *optional* [default = 0.5]
   	output grid resolution
    res : float, *optional* [default = 0.05]
   	spatial grid resolution to create a grid for the boundary pruning
    mxpts : float, *optional* [default = 1024]
   	maximum number of points allowed in a window
    minpts : float, *optional* [default = 16]
   	minimum number of points allowed in a window
    prc_overlap : float, *optional"  [default = 0]
        percentage overlap between windows

    Returns
    ----------
    self.data: list
   	list of M ndarrays, each containing n indices
        of original point cloud, toproc, to partition space
        to create M windows

    '''


### Detrend
    '''
    Detrend a Nx3 point cloud by specified method

    Syntax
    ----------
    detrended_pts = pysesa.detrend(points, proctype, res, method).getdata()

    Parameters
    ----------
    points : ndarray
   	Nx3 point cloud
    proctype : int
   	type of detrending.
        1 = remove mean
        2 = remove Ordinary least squares plane
        3 = remove Robust linear model plane
        4 = remove Orthogonal Distance Regression plane
        5 = remove Savitsky-Golay digital filter, order 1

    Other Parameters
    ----------
    res : float, *optional* [default = 0.05]
   	for proctype==4 only
        spatial grid resolution to create a grid
    method : str, *optional* [default = 'nearest']
   	for proctype==4 only
   	gridding type

    Returns
    ----------
    self.data: ndarray
   	Nx3 detrended point cloud

    '''

### Spatial
    '''
    Calculate spatial statistics of a Nx3 point cloud

    Syntax
    ----------
    stats = pysesa.spatial(points).getdata()
    centroids = pysesa.spatial(points).getcentroid()
    stats = pysesa.spatial(points).getstats()

    Parameters
    ----------
    points : ndarray
   	Nx3 point cloud

    Returns [requested through .getdata()]
    ----------
    self.data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
        z_mean = centroid in amplitude
        z_max = max amplitude
        z_min = min amplitude
        z_range = range in amplitude
        sigma = standard deviation of amplitudes
        skewness = skewness of amplitudes
        kurtosis = skewness of amplitudes
        n = number of 3D coordinates

    Returns [requested through .getcentroid()]
    ----------
    self.centroid: list
   	1x3 point cloud centroid [x,y,z]

    Returns [requested through .getstats()]
    ----------
    self.stats: list
        z_mean = centroid in amplitude
        z_max = max amplitude
        z_min = min amplitude
        z_range = range in amplitude
        sigma = standard deviation of amplitudes
        skewness = skewness of amplitudes
        kurtosis = skewness of amplitudes
        n = number of 3D coordinates

    '''


### Lengthscale
    '''
    Calculates the integral lengthscale of a Nx3 point cloud
    using 1 of 3 available methods
    and also returns the tapered 2D grid of 3D pointcloud for spectral analysis

    Syntax
    ----------
    im = pysesa.lengthscale(points, res, lentype, taper, method).getdata()
    lengthscale = pysesa.lengthscale(points, res, lentype, taper, method).getlengthscale()

    Parameters
    ----------
    points : ndarray
   	Nx3 point cloud

    Other Parameters
    ----------
    res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
    lentype : int, *optional* [default = 0, l<0.5]
   	lengthscale type:
        1, l<0.5
        2, l<1/e
        3, l<0
    taper : int, *optional* [default = Hanning]
   	flag for taper type:
        1, Hanning (Hann)
        2, Hamming
        3, Blackman
        4, Bartlett
    method : str, *optional* [default = 'nearest']
   	gridding type

    Returns [requested through .getdata()]
    ----------
    self.data: ndarray
   	tapered 2D grid of 3D pointcloud

    Returns [requested through .getlengthscale()]
    ----------
    self.lengthscale: float
   	integral lengthscale

    '''

### Spectral
    '''
    Calculate spectral statistics of a Nx3 point cloud

    Syntax
    ----------
    data = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getdata()
    lengths = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengths()
    psdparams= pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getstats()
    lengthscale = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengthscale()
    moments = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getmoments()

    Parameters
    ----------
    points : ndarray
   	Nx3 point cloud

    Other Parameters
    ----------
    nbin : int, *optional* [default = 20]
        number of bins for power spectral binning
    res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
    proctype : int, *optional* [default = 1, no spectral smoothing]
   	proctype type:
        1, no spectral smoothing
    lentype : int, *optional* [default = 1, l<0.5]
   	lengthscale type:
        1, l<0.5
        2, l<1/e
        3, l<0
    taper : int, *optional* [default = Hanning]
   	flag for taper type:
        1, Hanning (Hann)
        2, Hamming
        3, Blackman
        4, Bartlett
    method : str, *optional* [default = 'nearest']
   	gridding type

    Returns [requested through .getdata()]
    ----------
    self.data: list
   	slope = slope of regression line through log-log 1D power spectral density
        intercept = intercept of regression line through log-log 1D power spectral density
        r_value = correlation of regression through log-log 1D power spectral density
        p_value = probability that slope of regression through log-log 1D power spectral density is not zero
        std_err = standard error of regression through log-log 1D power spectral density
        d = fractal dimension
        l = integral lengthscale
        wmax = peak wavelength
        wmean = mean wavelength
        rms1 = RMS amplitude from power spectral density
        rms2 = RMS amplitude from bin averaged power spectral density
        Z = zero-crossings per unit length
        E = extreme per unit length
        sigma = RMS amplitude
        T0_1 = average spatial period (m_0/m_1)
        T0_2 = average spatial period (m_0/m_2)^0.5
        sw1 = spectral width
        sw2 = spectral width (normalised radius of gyration)
        m0 = zeroth moment of spectrum
        m1 = first moment of spectrum
        m2 = second moment of spectrum
        m3 = third moment of spectrum
        m4 = fourth moment of spectrum
        phi = effective slope (degrees)

    Returns [requested through .getpsdparams()]
    ----------
    self.psdparams: list
   	slope = slope of regression line through log-log 1D power spectral density
        intercept = intercept of regression line through log-log 1D power spectral density
        r_value = correlation of regression through log-log 1D power spectral density
        p_value = probability that slope of regression through log-log 1D power spectral density is not zero
        std_err = standard error of regression through log-log 1D power spectral density
        d = fractal dimension

    Returns [requested through .getlengths()]
    ----------
    self.lengths: list
        wmax = peak wavelength
        wmean = mean wavelength
        rms1 = RMS amplitude from power spectral density
        rms2 = RMS amplitude from bin averaged power spectral density

    Returns [requested through .getlengthscale()]
    ----------
    self.lengthscale: float
        l = integral lengthscale

    Returns [requested through .getmoments()]
    ----------
    self.moments: list
        Z = zero-crossings per unit length
        E = extreme per unit length
        sigma = RMS amplitude
        T0_1 = average spatial period (m_0/m_1)
        T0_2 = average spatial period (m_0/m_2)^0.5
        sw1 = spectral width
        sw2 = spectral width (normalised radius of gyration)
        m0 = zeroth moment of spectrum
        m1 = first moment of spectrum
        m2 = second moment of spectrum
        m3 = third moment of spectrum
        m4 = fourth moment of spectrum
        phi = effective slope (degrees)

    '''

### Process
    '''
    Calculate spectral and spatial statistics of a Nx3 point cloud

    Syntax
    ----------
    () = pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)

    Parameters
    ----------
    infile : str
   	ASCII file containing an Nx3 point cloud in 3 columns

    Other Parameters
    ----------
    out : float, *optional* [default = 0.5]
   	output grid resolution
    detrend : int, *optional* [default = 4]
   	type of detrending.
        1 = remove mean
        2 = remove Ordinary least squares plane
        3 = remove Robust linear model plane
        4 = remove Orthogonal Distance Regression plane
        5 = remove Savitsky-Golay digital filter, order 1
    proctype : int, *optional* [default = 1, no spectral smoothing]
   	proctype type:
        1 = spectral only, no spectral smoothing
        2 = spectral only, spectrum smoothed with Gaussian
        3 = spatial only
        4 = spatial + spectrum, no spectral smoothing
        5 = spatial + spectrum smoothed with Gaussian
    mxpts : float, *optional* [default = 1024]
   	maximum number of points allowed in a window
    res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
    nbin : int, *optional* [default = 20]
        number of bins for power spectral binning
    lentype : int, *optional* [default = 1, l<0.5]
   	lengthscale type:
        1 = l<0.5
        2 = l<1/e
        3 = l<0
    minpts : float, *optional* [default = 16]
   	minimum number of points allowed in a window
    taper : int, *optional* [default = Hanning]
   	flag for taper type:
        1 = Hanning (Hann)
        2 = Hamming
        3 = Blackman
        4 = Bartlett
    prc_overlap : float, *optional"  [default = 0]
        percentage overlap between windows


    Returns [proctype = 1 or proctype = 2]
    ----------
    data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
   	slope = slope of regression line through log-log 1D power spectral density
        intercept = intercept of regression line through log-log 1D power spectral density
        r_value = correlation of regression through log-log 1D power spectral density
        p_value = probability that slope of regression through log-log 1D power spectral density is not zero
        std_err = standard error of regression through log-log 1D power spectral density
        d = fractal dimension
        l = integral lengthscale
        wmax = peak wavelength
        wmean = mean wavelength
        rms1 = RMS amplitude from power spectral density
        rms2 = RMS amplitude from bin averaged power spectral density
        Z = zero-crossings per unit length
        E = extreme per unit length
        sigma = RMS amplitude
        T0_1 = average spatial period (m_0/m_1)
        T0_2 = average spatial period (m_0/m_2)^0.5
        sw1 = spectral width
        sw2 = spectral width (normalised radius of gyration)
        m0 = zeroth moment of spectrum
        m1 = first moment of spectrum
        m2 = second moment of spectrum
        m3 = third moment of spectrum
        m4 = fourth moment of spectrum
        phi = effective slope (degrees)

    Returns [proctype = 3]
    ----------
    data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
        z_mean = centroid in amplitude
        z_max = max amplitude
        z_min = min amplitude
        z_range = range in amplitude
        sigma = standard deviation of amplitudes
        skewness = skewness of amplitudes
        kurtosis = skewness of amplitudes
        n = number of 3D coordinates

    Returns [proctype = 4 or proctype = 4]
    ----------
    data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
        z_mean = centroid in amplitude
        z_max = max amplitude
        z_min = min amplitude
        z_range = range in amplitude
        sigma = standard deviation of amplitudes
        skewness = skewness of amplitudes
        kurtosis = skewness of amplitudes
        n = number of 3D coordinates
   	slope = slope of regression line through log-log 1D power spectral density
        intercept = intercept of regression line through log-log 1D power spectral density
        r_value = correlation of regression through log-log 1D power spectral density
        p_value = probability that slope of regression through log-log 1D power spectral density is not zero
        std_err = standard error of regression through log-log 1D power spectral density
        d = fractal dimension
        l = integral lengthscale
        wmax = peak wavelength
        wmean = mean wavelength
        rms1 = RMS amplitude from power spectral density
        rms2 = RMS amplitude from bin averaged power spectral density
        Z = zero-crossings per unit length
        E = extreme per unit length
        sigma = RMS amplitude
        T0_1 = average spatial period (m_0/m_1)
        T0_2 = average spatial period (m_0/m_2)^0.5
        sw1 = spectral width
        sw2 = spectral width (normalised radius of gyration)
        m0 = zeroth moment of spectrum
        m1 = first moment of spectrum
        m2 = second moment of spectrum
        m3 = third moment of spectrum
        m4 = fourth moment of spectrum
        phi = effective slope (degrees)
    '''

### Write
    '''
    Custom fast numpy array to comma-delimited ASCII txt file

    Syntax
    ----------
    () = pysesa.write.txtwrite(infile, towrite)

    Parameters
    ----------
    outfile : str
   	name of file to write to
    towrite : ndarray
   	ndarray containing Nx3 point cloud

    Other Parameters
    ----------
    header : str, *optional* [default = None]
   	header string

    Returns
    ----------
    None

    '''
