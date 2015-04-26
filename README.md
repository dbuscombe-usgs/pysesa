![alt tag](http://dbuscombe-usgs.github.io/figs/pysesa_colour.jpg)

PySESA - a Python framework for Spatially Explicit Spectral Analysis

PySESA is an open-source project dedicated to provide a generic Python framework 
for spatially explicit statistical analyses of point clouds and other geospatial data, 
in the spatial and frequency domains, for use in the geosciences

The program is detailed in:
Buscombe, D. "Spatially explicit spectral analysis of point clouds and geospatial data", forthcoming.

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
        | dbuscombe@usgs.gov
Version: 0.0.14    |  Revision: Apr, 2015

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

### Automatic Installation from PyPI 
The PyPI repository is here: https://pypi.python.org/pypi/pysesa

```
pip uninstall pysesa (removes previous installation)
pip install pysesa
```

Automatic Installation from github:

```
git clone git@github.com:dbuscombe-usgs/pysesa.git
cd pysesa
python setup.py install
```

or a local installation:

```
python setup.py install --user
```

or with admin privileges, e.g.:

```
sudo python setup.py install
```


### Notes for Linux users

You could try before you install, using a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install numpy
pip install cython
pip install scipy
pip install joblib
pip install statsmodels
pip install matplotlib
pip install ift_nifty
pip install pysesa
python -c "import pysesa; pysesa.test()"
deactivate (or source venv/bin/deactivate)
```

The results will live in "venv/lib/python2.7/site-packages/pysesa"


###Manual Installation

PYTHON LIBRARIES YOU MAY NEED TO INSTALL TO USE pysesa:

1. Nifty (http://www.mpa-garching.mpg.de/ift/nifty/index.html)
2. SciPy (http://www.scipy.org/scipylib/download.html)
3. Numpy (http://www.scipy.org/scipylib/download.html)
4. Matplotlib (http://matplotlib.org/downloads.html)
5. cython (http://cython.org/)
6. statsmodels (http://statsmodels.sourceforge.net/)
7. mayavi, used in some pysesa::plot modules (http://docs.enthought.com/mayavi/mayavi/)

All of the above are available through pip (https://pypi.python.org/pypi/pip) and easy_install (https://pythonhosted.org/setuptools/easy_install.html)

###Installation on Amazon Linux EC-2 instance
It's best to install numpy, scipy, cython and matplotlib through the OS package manager:

```
  sudo yum install gcc gcc-c++
  sudo yum install python27-numpy python27-Cython python27-scipy python27-matplotlib
```   

Then pysesa using pip (which will install nifty, joblib and statsmodels):

```
  sudo pip install pysesa
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

   out = 1 #m output grid
   detrend = 4 #ODR plane
   proctype = 1 #Processing spectral parameters (no smoothing)
   mxpts = 1024 # max pts per window
   res = 0.05 #cm internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # Hann taper
   prc_overlap = 50 # no overlap between successive windows
   minpts = 64 # min pts per window

   pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)

``` 

### Support

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 


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
12. plot: some plotting utilities for outputs in 2d and 3d

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

### Sgolay
    '''
    Create a Savitsky-Golay digital filter from a 2D signal
    based on code from http://www.scipy.org/Cookbook/SavitzkyGolay

    Syntax
    ----------
    Z = pysesa.sgolay(z, window_size, order).getdata()

    Parameters
    ----------
    z : array_like, shape (N,)
      the 2D signal.
    window_size : int
       the length of the window. Must be an odd integer number.
    order : int
       the order of the polynomial used in the filtering.
       Must be less than `window_size` - 1.

    Returns
    -------
    self.data : ndarray, shape (N)
       the smoothed signal.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
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
        2, spectrum smoothed with Gaussian
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


### Plot
    '''
    Initialise the pysesa plot class
    Takes a file output from pysesa.process and allows a number
    of different 2d or 3d plots of the outputs

    Syntax
    ----------
    p = pysesa.plot()
    p = pysesa.plot('/home/my_pysesa_output_file.xyz')

    Parameters
    -----------
    pysesa_file : str
   	pysesa::process output file

    If no arguments given, it prompts you to choose a pysesa::process output file

    Returns
    ----------
    self : instance
      pysesa.plot instance

    DATA functions
    ---------------
    pc = p.get_pc() : ndarray
      NxM contents of pysesa_file
  
    xyz = p.get_xyz() : ndarray        
      Nx3 contents of raw point cloud (the file processed by pysesa::process)

    vars = p.parse_pc_vars() : dict
      NxM contents of pysesa_file parsed into dict object
      1 key per variable in p.get_pc()
      vars.keys() returns list of variables in dict

    2D plotting functions (mayavi not required)
    --------------------------------------------

    p.grd_xyz()
    ------------- 
      produces 2d plot of the gridded [x,y,z] surface made from decimated point cloud,
      as returned by parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_xyz(res, azimuth, altitude, zf, cmap, dpi, alpha, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

    p.grd_var()
    ------------      
      produces 2d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()
      e.g. p.grd_var('sigma')

      Syntax
      ----------
      [] = p.grd_var(var, res, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var('sigma')


      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

      filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

      alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels


    p.grd_vars()
    -------------
      produces a 2d plot of the gridded surface made from each output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_vars(res, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

      filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

      alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels


    3D plotting functions (mayavi not required)
    --------------------------------------------

    p.plt_xyz()
    ------------ 
      produces 3d plot of Nx3 contents of raw point cloud, as returned by p.get_xyz()

      Syntax
      ----------
      [] = p.plt_xyz(elev, azim, markersize, dpi, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      elev : float, *optional* [default = 65]
           the elevation angle in the z plane

      azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane 

      markersize : float, *optional* [default = 0.01]
           marker size in x and y axes units

      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels       


    p.plt_xy_var()  
    --------------- 
      produces 3d plot of 1 output variable in p.parse_pc_vars(), e.g. p.grd_var('sigma')

      Syntax
      ----------
      [] = p.plt_xy_var(var, log_scale, dpi, markersize, ticksize, labelsize, elev, azim)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var('sigma')

      Optional Parameters
      --------------------
      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      markersize : float, *optional* [default = 5]
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

      elev : float, *optional* [default = 65]
           the elevation angle in the z plane

      azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane          


    p.plt_xy_vars()
    ----------------
      produces a 3d plot of each output variable in p.parse_pc_vars() 

      Syntax
      ----------
      [] = p.plt_xy_vars(log_scale, dpi, markersize, ticksize, labelsize, elev, azim)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      markersize : float, *optional* [default = 5]
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

      elev : float, *optional* [default = 65]
           the elevation angle in the z plane

      azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane        
   

    3D plotting functions (requires mayavi) 
    ----------------------------------------

    p.grd_xyz3d() 
    -------------  
      produces 3d plot of the gridded surface made from the Nx3 contents of raw point cloud,
      as returned by p.get_xyz()

      Syntax
      ----------
      [] = p.grd_xyz3d(res, cmap, pitch, azimuth, distance, xsize, ysize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()


      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

      ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension
   
    p.grd_var_3d()  
    ---------------
      produces 3d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_var_3d(var, res, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var_3d('sigma')


      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

      filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

      xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

      ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension
          
    p.grd_vars_3d() 
    ----------------
      produces a 3d plot of the gridded surface made from each output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_vars_3d(res, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

      smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

      filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

      xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

      ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension

    '''
