.. _examples:

Examples
===================

Process module
---------------
The submodule *process* allows full control over all types of workflows through use of a number of processing flags. A minimum working example usage of the the PySESA module, accepting all default values for parameters, is::

  import pysesa
  infile = `/home/me/mypointcloudfile.txt'
  pysesa.process(infile)

This instance writes out the following results file whose name contains some of the processing parameters::
 
  /home/me/mypointcloudfile.txt_zstat_detrend4_outres0.5_proctype1_mxpts512_minpts16.xyz

Test module
---------------
The above is the same as passing a list of default-valued variables to *process*, which is included for completeness in the *test* module::

  out = 1 		# 1 m output grid
  detrend = 4 		# detrend type: ODR plane
  proctype = 1 		# Processing type: spectral parameters (no smoothing) only
  mxpts = 1024 		# Maximum points per window
  res = 0.05 		# 5 cm grid resolution for detrending and spectral analysis
  nbin = 20 		# Number of bins for spectral binning
  lentype = 1 		# Integral lengthscale type: l<0.5
  taper = 1 		# Hann taper before spectral analysis
  prc_overlap = 0 	# No overlap between successive windows
  minpts = 64 		# Minimum points per window

  pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)

Minimal example on 1 window
----------------------------
A minimal example analysis of spatial and spectral analysis on just 1 window of data::

  # import module
  import pysesa

  # read point cloud from file
  pointcloud = pysesa.read.txtread(infile)

  # create windows of data
  windows = pysesa.partition(pointcloud).getdata()

  # process window number 50
  k=50

  # get all spectral statistics for that window
  spec_stats = pysesa.spectral(pointcloud[windows[k],:3].astype('float64')).getdata()

  # get all spatial statistics for that window
  spat_stats = pysesa.spatial(pointcloud[windows[k],:3].astype('float64')).getdata()


Minimal example, all windows using parallel processing
--------------------------------------------------------
and to extend this to all windows, utilising parallel processing over all available cores, could be achieved using the following minimal example::
 
  # define a function that will get repeatedly read by the parallel processing queue
  def get_spat_n_spec(pts):
     return pysesa.spatial(pts.astype('float64')).getdata() + pysesa.spectral(pts.astype('float64')).getdata()

  # import the parallel processing libraries
  from joblib import Parallel, delayed, cpu_count

  # Processing type: spatial plus spectral parameters (no smoothing)
  proctype = 4

  # process each window with all available cores, by queueing each window in a sequence
  # and processing until they are all done
  w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spat_n_spec)(pointcloud[windows[k],:3])
	for k in xrange(len(windows)))

  # parse out the outputs into variables
  x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, ...
  slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, ...
  Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi = zip(*w)

Lengthscale module
-------------------
To obtain just the integral lengthscale of the kth window, detrended using the orthogonal distance regression detrending technique, one could use::

  detrend = 4 	# Orthogonal distance regression
  pysesa.lengthscale(pysesa.detrend(pointcloud[windows[k],:3],detrend).getdata()).getlengthscale()

Spatial module
---------------
and to get the spatial statistics from the same data::

  pysesa.spatial(pysesa.detrend(pointcloud[windows[k],:3],detrend).getdata()).getdata()

Spectral module
------------------
Here, the output grid resolution is changed to 25 cm, and the various outputs from the *spectral* module are obtained separately::

  # 25 cm output grid
  out = 0.25

  # re-create windows of data
  windows = pysesa.partition(pointcloud, out).getdata()

  result = pysesa.spectral(pointcloud[windows[k],:3].astype('float64'))

  # get all spectral parameters
  result.getdata()

  # get the fit parameters for log-log power spectrum
  result.getpsdparams()

  # get integral lengthscale
  result.getlengthscale()

  # get spectral moment parameters
  result.getmoments()

  # get rms and wavelength parameters
  result.getlengths()

Plot module
------------------
This assumes you have run the *process* module and have an output file ('/home/my_pysesa_output_file.xyz')::

  # load pysesa
  import pysesa

  # create a pysesa::plot instance
  p = pysesa.plot('/home/my_pysesa_output_file.xyz')

  # create a 3d plot of the point cloud
  p.plt_xyz()

  # create a 2d plot of the gridded surface from the decimated point cloud
  p.grd_xyz()

  # create a 3d plot of the gridded surface from the decimated point cloud
  # colour-coded by amplitude
  p.grd_xyz3d()

  # create a 3d plot of the decimated spectral slope
  # colour-coded by amplitude
  p.plt_xy_var('slope')

  # create a 3d plot of all output decimated parameters
  # colour-coded by amplitude
  p.plt_xy_vars()

  # create a 3d plot of decimated fractal dimension
  # gridded and colour-coded by amplitude
  p.grd_var_3d('d')

  # create a 3d plot of all output decimated parameters
  # gridded and colour-coded by amplitude
  p.grd_vars_3d()

  # plot also supports data retrieval
  # retrieve the original point cloud
  xyz = p.get_xyz()

  # retrieve the decimated point cloud of all parameters
  pc = p.get_pc()

  # retrieve the decimated point cloud of all parameters
  # in dictionary format
  pc_dict = p.parse_pc_vars()
  # then show what's in there
  print pc.dict.keys()

.. image:: _static/pysesa_colour.jpg


