## PySESA (Python program for Spatially Explicit Spectral Analysis)
## has been developed at the Grand Canyon Monitoring & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/pysesa>
##
##This software is in the public domain because it contains materials that originally came from
##the United States Geological Survey, an agency of the United States Department of Interior.
##For more information, see the official USGS copyright policy at
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
PySESA - a Python framework for Spatially Explicit Spectral Analysis

PySESA is an open-source project dedicated to provide a generic Python framework
for spatially explicit statistical analyses of point clouds and other geospatial data,
in the spatial and frequency domains, for use in the geosciences

The program is detailed in:
   Buscombe, D. (2016) "Computational considerations for spatially explicit spectral analysis of point clouds and geospatial data", 86, 92-108, 10.1016/j.cageo.2015.10.004.

:Author:
    Daniel Buscombe
    Grand Canyon Monitoring and Research Center
    United States Geological Survey
    Flagstaff, AZ 86001
    dbuscombe@usgs.gov

For more information visit http://dbuscombe-usgs.github.io/pysesa/

:install:
    python setup.py install
    sudo python setup.py install

:test:
    python -c "import pysesa; pysesa.test()"

:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information,
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only
    and does not imply endorsement by the U.S. government.

"""

# =========================================================
# ======================== libraries ======================
# =========================================================

from __future__ import division
import numpy as np
from joblib import Parallel, delayed, cpu_count
from time import clock, time
import os #, sys, getopt

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

import pysesa

import warnings
warnings.filterwarnings("ignore")

# =========================================================
# ===================== subfunctions ======================
# =========================================================

# =========================================================
def ascol( arr ):
   '''
   reshapes row matrix to be a column matrix (N,1).
   '''
   if len( arr.shape ) == 1: arr = arr.reshape( ( arr.shape[0], 1 ) )
   return arr

#==================================================
def get_spec(pts, spectype, out, detrend, res, method, nbin, lentype, taper):
   '''
   call the spectral analysis routine for detrended window of point cloud
   Gets called by the parallel processing queue
   '''
   try:
      pts = pysesa.detrend(pts, detrend, res, method).getdata()
      return pysesa.spatial(pts).getcentroid() + pysesa.spectral(pts, nbin, res, spectype, lentype, taper, method).getdata()
   except:
      return [np.ones(27)*np.nan]

#==================================================
def get_spat(pts, detrend, res, method):
   '''
   call the spatial analysis routine for detrended window of point cloud
   Gets called by the parallel processing queue
   '''
   try:
      pts = pysesa.detrend(pts, detrend, res, method).getdata()
      return pysesa.spatial(pts).getdata()
   except:
      return [np.ones(10)*np.nan]

#==================================================
def get_spec_spat(pts, spectype, out, detrend, res, method, nbin, lentype, taper):
   '''
   call the spectral and spatial analysis routine for detrended window of point cloud
   Gets called by the parallel processing queue
   '''
   try:
      pts = pysesa.detrend(pts, detrend, res, method).getdata()
      return pysesa.spatial(pts).getdata() + pysesa.spectral(pts, nbin, res, spectype, lentype, taper, method).getdata()
   except:
      return [np.ones(34)*np.nan]

# =========================================================
# ==================== begin program ======================
# =========================================================
def process(infile, out=1, detrend=4, proctype=1, mxpts=1024, res=0.05, nbin=20, lentype=1, minpts=64, taper=1, prc_overlap=0, nchunks=1, filt=0): # bp=0
   '''
   Calculate spectral and spatial statistics of a Nx3 point cloud

   Syntax
   ----------
   () = pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap, nchunks, filt)

   Parameters
   -----------
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
   nchunks : int, *optional"  [default = 1]
        split data into nchunks and process each separately
        use only if receiving memory errors with very large datasets
   filt : int, *optional"  [default = 0]
        if filt==1, point cloud will be filtered prior to analysis
        using a simple thresholded standard deviation approach

   Returns [proctype = 1 or proctype = 2]
   ---------------------------------------
   data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
        z = centroid in vertical coordinate
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
   -------------------------
   data: list
   	x = centroid in horizontal coordinate
        y = centroid in laterial coordinate
        z = centroid in vertical coordinate
        z_mean = centroid in amplitude
        z_max = max amplitude
        z_min = min amplitude
        z_range = range in amplitude
        sigma = standard deviation of amplitudes
        skewness = skewness of amplitudes
        kurtosis = skewness of amplitudes
        n = number of 3D coordinates

   Returns [proctype = 4 or proctype = 5]
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
        phi = effective slope (degrees)

   '''

   print """
       ..  _ __  _   _   ___  ___  ___  __ _
       .. | '_ \| | | | / __|/ _ \/ __|/ _` |
       .. | |_) | |_| | \__ \  __/\__ \ (_| |
       .. | .__/ \__, | |___/\___||___/\__,_|
       .. |_|    |___/
       ..
       .. +-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
       .. |b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
       .. +-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
       .. +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
       .. |d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
       .. +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
       .. +-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
       .. |U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
       .. +-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
   """

   if infile:
      print 'Input file is %s' % (infile)

   if prc_overlap:
      prc_overlap = np.asarray(prc_overlap,int)
      print 'Percent overlap is %s' % (str(prc_overlap))
   elif prc_overlap==0:
      prc_overlap = np.asarray(prc_overlap,int)
      print 'Percent overlap is %s' % (str(prc_overlap))

   if out:
      out = np.asarray(out,float)
      print 'Output grid size is %s' % (str(out))

   if detrend:
      detrend = np.asarray(detrend,int)

      if detrend==1:
         print 'Detrend type: remove mean'
      if detrend==2:
         print 'Detrend type: Ordinary least squares plane'
      if detrend==3:
         print 'Detrend type: Robust linear model plane'
      if detrend==4:
         print 'Detrend type: Orthogonal Distance Regression plane'
      if detrend==5:
         print 'Detrend type: Savitsky-Golay digital filter, order 1'

   if proctype:
      proctype = np.asarray(proctype,int)
      if proctype==0:
         proctype = 1

      if proctype==1:
         print 'spectral parameters (no smoothing)'
      elif proctype==2:
         print 'spectral parameters (with smoothing)'
      elif proctype==3:
         print 'spatial parameters'
      elif proctype==4:
         print 'spatial parameters + spectral parameters (no smoothing)'
      elif proctype==5:
         print 'spatial parameters + spectral parameters (smoothing)'

   if res:
      res = np.asarray(res,float)
      print 'Res. is %s' % (str(res))

   if mxpts:
      mxpts = np.asarray(mxpts,int)
      print 'Max points per window is %s' % (str(mxpts))

   if minpts:
      minpts = np.asarray(minpts,int)
      print 'Min points per window is %s' % (str(minpts))

   if nbin:
      nbin = np.asarray(nbin,float)
      print 'Number of bins is %s' % (str(nbin))

   if lentype:
      lentype = int(lentype)
      if lentype==1:
         print "lengthscale type: l<0.5"
      elif lentype==2:
         print "lengthscale type: l<1/e"
      elif lentype==3:
         print "lengthscale type: l<0"

   if taper:
      taper = np.asarray(taper,int)
      if taper==1:
         print 'Hanning taper'
      elif taper==2:
         print 'Hamming taper'
      elif taper==3:
         print 'Blackman taper'
      else:
         print 'Bartlett taper'

   if nchunks:
      nchunks = np.asarray(nchunks,int)
      print 'Number of chunks to process separately is %s' % (str(nchunks))

   if filt:
      filt = np.asarray(filt,int)
      if filt==1:
         print 'Point cloud will be filtered'

   # start timer
   if os.name=='posix': # true if linux/mac or cygwin on windows
      start1 = time()
   else: # windows
      start1 = clock()

   method = 'nearest'

   #internal parameters for filtering
   k = 10 #number of neighbours
   std_dev = 2 #standard deviation multiplier
   n_iter = 2 #number of iterations

   #==============================================================================
   print "(1) Reading data from file ..."
   
   # check first for laz/las format
   if 'las' in infile[-3:]:
      toproc_init = pysesa.read.lasread(infile)
   elif 'laz' in infile[-3:]:
      toproc_init = pysesa.read.lasread(infile)
   else: # read in ascii 3-column file containing point cloud data
      toproc_init = pysesa.read.txtread(infile)

   #==============================================================================
   # if requested, filter data
   if filt==1:
       print "(1b) Filtering data ..."
       print "Size of original data: %s" % (str(len(toproc_init)))
       # initial pass
       _, toproc_init_f = pysesa.filter.filt_stdev(toproc_init, k = k, std_dev = std_dev)
       del toproc_init
       #iterate through n_iter to refine filtering
       for nn in range(n_iter):
           _, toproc_init_f = pysesa.filter.filt_stdev(toproc_init_f, k = k, std_dev = std_dev)
       toproc_init = np.copy(toproc_init_f)
       del toproc_init_f
       print "Size of filtered data: %s" % (str(len(toproc_init)))

       infile = infile.split('.')[-2]+'_filt.xyz'
       print "Writing filtered data to file: "+infile
       with open(infile, 'wb') as f:
          np.savetxt(f, toproc_init[np.where(toproc_init[:,-1])[0],:], fmt=' '.join(['%8.6f,'] * np.shape(toproc_init)[1])[:-1])

   #==============================================================================
   toproc2 = np.array_split(toproc_init, nchunks)
   del toproc_init

   ## number of points, undecimated
   orig_pts = len(np.vstack(toproc2))

   TOWRITE = []

   counter = 1
   for toproc in toproc2:

      print "Working on chunk %s out of %s chunks ... " % (str(counter), str(len(toproc2)))
      counter += 1

      #==============================================================================
      print "(2) Partitioning data into windows ... "
      # get indices to windows
      nr_pts = pysesa.partition(toproc, out, mxpts, minpts, prc_overlap).getdata() #res, bp

      #==============================================================================
      print "(3) Processing in parallel using %s processors ... " % (str(cpu_count()))

      #==============================================================================
      if (proctype==1) or (proctype==2):

		  #spectral, no smooth
		  if proctype==1:
			 try: #parallel processing with all available cores
				w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spec)(toproc[nr_pts[k],:3], 1, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))
			 except: #fall back to serial
				w = Parallel(n_jobs=1, verbose=0)(delayed(get_spec)(toproc[nr_pts[k],:3], 1, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))

		  #spectral, with smooth
		  if proctype==2:
			 try: #parallel processing with all available cores
				w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spec)(toproc[nr_pts[k],:3], 2, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))
			 except: #fall back to serial
				w = Parallel(n_jobs=1, verbose=0)(delayed(get_spec)(toproc[nr_pts[k],:3], 2, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))

		  try:
			 x, y, z, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi = zip(*w)
		  except:
			 w2 = []
			 for k in xrange(len(w)):
				if len(w[k])==27:
				   w2.append(w[k])
			 x, y, z, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi = zip(*w2)
			 del w2

		  del w

		  # combine into single matrix for writing to file
		  towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)),ascol(np.asarray(z)),ascol(np.asarray(slope)),ascol(np.asarray(intercept)),ascol(np.asarray(r_value)),ascol(np.asarray(p_value)),ascol(np.asarray(std_err)),ascol(np.asarray(d)),ascol(np.asarray(l)),ascol(np.asarray(wmax)),ascol(np.asarray(wmean)),ascol(np.asarray(rms1)),ascol(np.asarray(rms2)),ascol(np.asarray(Z)),ascol(np.asarray(E)),ascol(np.asarray(sigma)),ascol(np.asarray(T0_1)),ascol(np.asarray(T0_2)),ascol(np.asarray(sw1)),ascol(np.asarray(sw2)),ascol(np.asarray(m0)),ascol(np.asarray(m1)),ascol(np.asarray(m2)),ascol(np.asarray(m3)),ascol(np.asarray(m4)),ascol(np.asarray(phi)) ))

		  del x, y, z, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi

		  # remove rows with any NaNs
		  towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

		  TOWRITE.append(towrite)
		  del towrite

		  # make a header string for the output file
		  header = 'x, y, z, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi'

      #==============================================================================
      elif proctype==3: #spatial only

		  try: #parallel processing with all available cores
			 w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spat)(toproc[nr_pts[k],:3], detrend, res, method) for k in xrange(len(nr_pts)))
		  except: #fall back to serial
			 w = Parallel(n_jobs=1, verbose=0)(delayed(get_spat)(toproc[nr_pts[k],:3], detrend, res, method) for k in xrange(len(nr_pts)))

		  try:
			 x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n = zip(*w)
		  except:
			 w2 = []
			 for k in xrange(len(w)):
				if len(w[k])==10:
				   w2.append(w[k])
			 x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n = zip(*w2)
			 del w2

		  del w

		  # combine into single matrix for writing to file
		  towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)),ascol(np.asarray(z_mean)),ascol(np.asarray(z_max)),ascol(np.asarray(z_min)),ascol(np.asarray(z_range)),ascol(np.asarray(sigma)),ascol(np.asarray(skewness)),ascol(np.asarray(kurtosis)), ascol(np.asarray(n)) ))

		  del x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n

		  # remove rows with any NaNs
		  towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

		  TOWRITE.append(towrite)
		  del towrite

		  # make a header string for the output file
		  header = 'x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n'

      #==============================================================================
      elif (proctype==4) or (proctype==5): #spectral and spatial

		  #spatial + spectral, no smooth
		  if proctype==4:
			 try: #parallel processing with all available cores
				w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spec_spat)(toproc[nr_pts[k],:3], 1, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))
			 except: #fall back to serial
				w = Parallel(n_jobs=1, verbose=0)(delayed(get_spec_spat)(toproc[nr_pts[k],:3], 1, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))

		  #spatial + spectral, with smooth
		  if proctype==5:
			 try: #parallel processing with all available cores
				w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_spec_spat)(toproc[nr_pts[k],:3], 2, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))
			 except: #fall back to serial
				w = Parallel(n_jobs=1, verbose=0)(delayed(get_spec_spat)(toproc[nr_pts[k],:3], 2, out, detrend, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))

		  try:
			 x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi = zip(*w)
		  except:
			 w2 = []
			 for k in xrange(len(w)):
				if len(w[k])==34:
				   w2.append(w[k])
			 x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi = zip(*w2)
			 del w2

		  del w

		  # combine into single matrix for writing to file
		  towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)), ascol(np.asarray(z_mean)),ascol(np.asarray(z_max)),ascol(np.asarray(z_min)),ascol(np.asarray(z_range)),ascol(np.asarray(sigma)),ascol(np.asarray(skewness)),ascol(np.asarray(kurtosis)), ascol(np.asarray(n)), ascol(np.asarray(slope)),ascol(np.asarray(intercept)),ascol(np.asarray(r_value)),ascol(np.asarray(p_value)),ascol(np.asarray(std_err)),ascol(np.asarray(d)),ascol(np.asarray(l)),ascol(np.asarray(wmax)),ascol(np.asarray(wmean)),ascol(np.asarray(rms1)),ascol(np.asarray(rms2)),ascol(np.asarray(Z)),ascol(np.asarray(E)),ascol(np.asarray(sigma)),ascol(np.asarray(T0_1)),ascol(np.asarray(T0_2)),ascol(np.asarray(sw1)),ascol(np.asarray(sw2)),ascol(np.asarray(m0)),ascol(np.asarray(m1)),ascol(np.asarray(m2)),ascol(np.asarray(m3)),ascol(np.asarray(m4)),ascol(np.asarray(phi)) ))

		  del x, y, z_mean, z_max, z_min, z_range,  skewness, kurtosis, n #,sigma
		  del slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi

		  # remove rows with any NaNs
		  towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

		  TOWRITE.append(towrite)
		  del towrite

		  # make a header string for the output file
		  header = 'x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4, phi'

   towrite = np.vstack(TOWRITE)

   #==============================================================================
   print "(4) Writing data to file ..."

   # create a string for the output file
   outfile = infile+'_zstat_detrend'+str(detrend)+'_outres'+str(out)+'_proctype'+str(proctype)+'_mxpts'+str(mxpts)+'_minpts'+str(minpts)+'.xyz'

   try:
      # write the data to the file
      pysesa.write.txtwrite(outfile, towrite, header)

   except:
      with open(outfile, 'wb') as f:
         np.savetxt(f, towrite[np.where(towrite[:,-1])[0],:], header = header, fmt=' '.join(['%8.6f,'] * np.shape(towrite)[1])[:-1])

   x = np.copy(towrite)[:,6:]
   x_normed = (x - x.min(0)) / x.ptp(0)
   towrite2 = np.hstack((towrite[:,:6], x_normed))
   outfile = infile+'_zstat_detrend'+str(detrend)+'_outres'+str(out)+'_proctype'+str(proctype)+'_mxpts'+str(mxpts)+'_minpts'+str(minpts)+'_norm1.xyz'

   try:
      # write the data to the file
      pysesa.write.txtwrite(outfile, towrite2, header)

   except:
      with open(outfile, 'wb') as f:
         np.savetxt(f, towrite2[np.where(towrite2[:,-1])[0],:], header = header, fmt=' '.join(['%8.6f,'] * np.shape(towrite2)[1])[:-1])

   towrite2 = np.hstack((towrite[:,:6], x_normed*255))
   outfile = infile+'_zstat_detrend'+str(detrend)+'_outres'+str(out)+'_proctype'+str(proctype)+'_mxpts'+str(mxpts)+'_minpts'+str(minpts)+'_norm255.xyz'

   try:
      # write the data to the file
      pysesa.write.txtwrite(outfile, towrite2, header)

   except:
      with open(outfile, 'wb') as f:
         np.savetxt(f, towrite2[np.where(towrite2[:,-1])[0],:], header = header, fmt=' '.join(['%8.6f,'] * np.shape(towrite2)[1])[:-1])

   # stop the clock
   if os.name=='posix': # true if linux/mac
      elapsed = (time() - start1)
   else: # windows
      elapsed = (clock() - start1)

   print "Done! %s points decimated to %s points. Program ran for %s seconds" % (str(orig_pts), str(len(towrite)), str(elapsed))
