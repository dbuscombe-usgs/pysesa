.. _modules:

***************
Modules
***************

.. image:: _static/pysesa_colour.jpg

.. _overview:

Overview
=========

The programs in this package are as follows:

**read**: read a 3-column space, comma or tab delimited text file

**partition**: partition a Nx3 point cloud into M windows of nx3 points with specified spacing between centroids of adjacent windows and with specified overlap between windows.

**detrend**: returns detrended amplitudes of a Nx3 point cloud

**sgolay**: returns the Savitsky-Golay digital filter of a 2D signal

**spatial**: calculate spatial statistics of a Nx3 point cloud

**lengthscale**: calculates the integral lengthscale of a Nx3 point cloud

**spectral**: calculate spectral statistics of a Nx3 point cloud

**process**: allows control of inputs to all modules (full workflow)

**write**: write program outputs to a comma delimited text file 

**test**: program testing suite 

These are all command-line/modular programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options


.. _process:

Main program: process
=====================

Calculate spectral and spatial statistics of a Nx3 point cloud


Syntax
----------

You call the function like this::

  () = pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)


Parameters
------------
   infile : str
   	ASCII file containing an Nx3 point cloud in 3 columns

Other Parameters
-----------------
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
----------------------------------------
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


Returns [proctype = 3]
------------------------
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
-----------------------------------------
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


  .. image:: _static/pysesa_colour.jpg

