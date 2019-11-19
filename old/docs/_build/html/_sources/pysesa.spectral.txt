.. _pysesa.spectral:

pysesa.spectral module
=======================

Calculate spectral statistics of a Nx3 point cloud

Syntax
----------

You call the function like this::

  data = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getdata()

or::

  lengths = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengths()

or::

  psdparams= pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getstats()

or::

  lengthscale = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengthscale()

or::

  moments = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getmoments()

Parameters
------------
points : ndarray
   	Nx3 point cloud

Other Parameters
-------------------
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
----------------------------------------
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
--------------------------------------------
self.psdparams: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

Returns [requested through .getlengths()]
-------------------------------------------
self.lengths: list
        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density


Returns [requested through .getlengthscale()]
----------------------------------------------
   self.lengthscale: float
        l = integral lengthscale

Returns [requested through .getmoments()]
------------------------------------------
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

  .. image:: _static/pysesa_colour.jpg
