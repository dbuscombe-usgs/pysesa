## PySESA (Python program for Spatially Explicit Spectral Analysis) 
## has been developed at the Grand Canyon Monitorinf & Research Center,
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

# =========================================================
# ======================== libraries ======================
# =========================================================

from __future__ import division
import numpy as np
from joblib import Parallel, delayed, cpu_count
from time import clock, time
import os, sys, getopt

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

#import pysesa_lengthscale
import pysesa.spatial 
import pysesa.spectral 
import pysesa.partition 
import pysesa.read 
import pysesa.write 
import pysesa.detrend 

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
   pts = pysesa.detrend.detrend(pts, detrend, res, method).getdata()
   return pysesa.spatial.spatial(pts).getcentroid() + pysesa.spectral.spec(pts, nbin, res, spectype, lentype, taper, method).getdata()

#==================================================
def get_spat(pts, detrend, res, method):
   '''
   call the spatial analysis routine for detrended window of point cloud
   Gets called by the parallel processing queue
   '''
   pts = pysesa.detrend.detrend(pts, detrend, res, method).getdata()
   return pysesa.spatial.spatial(pts).getdata()

#==================================================
def get_spec_spat(pts, spectype, out, detrend, res, method, nbin, lentype, taper):
   '''
   call the spectral and spatial analysis routine for detrended window of point cloud
   Gets called by the parallel processing queue
   '''
   pts = pysesa.detrend.detrend(pts, detrend, res, method).getdata()
   return pysesa.spatial.spatial(pts).getdata() + pysesa.spectral.spec(pts, nbin, res, spectype, lentype, taper, method).getdata()

# =========================================================
# ==================== begin program ======================
# =========================================================
def process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap):

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

   if not infile:
      infile = 'example_100000pts.xyz'
      print '[Default] Input file is %s' % (infile)

   if prc_overlap:
      prc_overlap = np.asarray(prc_overlap,float)
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
      if lentype==0:
         print "lengthscale type: l<0"
      elif lentype==1:
         print "lengthscale type: l<1/e"
      else:
         print "lengthscale type: l<0.5"

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

   if not prc_overlap:
      prc_overlap = 0
      print '[Default] Percentage overlap is %s' % (str(prc_overlap))

   if not out:
      out = 0.5
      print '[Default] Output grid size is %s' % (str(out))

   if not detrend and detrend!=0:
      detrend = 4
      print '[Default] Type of detrend is %s (ODR plane)' % (str(detrend))

   if not proctype:
      proctype = 1
      print '[Default] Processing spectral parameters (no smoothing)'

   if not res:
      res = 0.05
      print '[Default] Res. is %s' % (str(res))

   if not mxpts:
      mxpts = 512
      print '[Default] Max points per window is %s' % (str(mxpts))

   if not minpts:
      minpts = 16
      print '[Default] Min points per window is %s' % (str(minpts))

   if not nbin:
      nbin = 20
      print '[Default] Number of bins is %s' % (str(nbin))

   if not lentype:
      lentype = 0
      print "[Default] lengthscale type: l<0.5"

   if not taper:
      taper = 1
      print '[Default] Hanning taper'

   # start timer
   if os.name=='posix': # true if linux/mac or cygwin on windows
      start1 = time()
   else: # windows
      start1 = clock()

   method = 'nearest'

   #==============================================================================
   print "(1) Reading data from file ..."
   # read in ascii 3-column file containing point cloud data
   toproc = pysesa.read.txtread(infile)

   ## number of points, undecimated
   orig_pts = len(toproc)

   #==============================================================================
   print "(2) Partitioning data into windows ... " 
   # get indices to windows
   nr_pts = pysesa.partition.partition(toproc, out, res, mxpts, minpts, prc_overlap).getdata()

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
         x, y, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = zip(*w)
      except:
         w2 = []
         for k in xrange(len(w)):
            if len(w[k])==25:
               w2.append(w[k])
         x, y, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = zip(*w2)
         del w2   

      del w

      # combine into single matrix for writing to file
      towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)),ascol(np.asarray(slope)),ascol(np.asarray(intercept)),ascol(np.asarray(r_value)),ascol(np.asarray(p_value)),ascol(np.asarray(std_err)),ascol(np.asarray(d)),ascol(np.asarray(l)),ascol(np.asarray(wmax)),ascol(np.asarray(wmean)),ascol(np.asarray(rms1)),ascol(np.asarray(rms2)),ascol(np.asarray(Z)),ascol(np.asarray(E)),ascol(np.asarray(sigma)),ascol(np.asarray(T0_1)),ascol(np.asarray(T0_2)),ascol(np.asarray(sw1)),ascol(np.asarray(sw2)),ascol(np.asarray(m0)),ascol(np.asarray(m1)),ascol(np.asarray(m2)),ascol(np.asarray(m3)),ascol(np.asarray(m4)) ))

      del x, y, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4

      # remove rows with any NaNs
      towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

      # make a header string for the output file
      header = 'x, y, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4'

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
         x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = zip(*w)
      except:
         w2 = []
         for k in xrange(len(w)):
            if len(w[k])==33:
               w2.append(w[k])
         x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = zip(*w2)
         del w2   

      del w

      # combine into single matrix for writing to file
      towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)), ascol(np.asarray(z_mean)),ascol(np.asarray(z_max)),ascol(np.asarray(z_min)),ascol(np.asarray(z_range)),ascol(np.asarray(sigma)),ascol(np.asarray(skewness)),ascol(np.asarray(kurtosis)), ascol(np.asarray(n)), ascol(np.asarray(slope)),ascol(np.asarray(intercept)),ascol(np.asarray(r_value)),ascol(np.asarray(p_value)),ascol(np.asarray(std_err)),ascol(np.asarray(d)),ascol(np.asarray(l)),ascol(np.asarray(wmax)),ascol(np.asarray(wmean)),ascol(np.asarray(rms1)),ascol(np.asarray(rms2)),ascol(np.asarray(Z)),ascol(np.asarray(E)),ascol(np.asarray(sigma)),ascol(np.asarray(T0_1)),ascol(np.asarray(T0_2)),ascol(np.asarray(sw1)),ascol(np.asarray(sw2)),ascol(np.asarray(m0)),ascol(np.asarray(m1)),ascol(np.asarray(m2)),ascol(np.asarray(m3)),ascol(np.asarray(m4)) ))

      del x, y, z_mean, z_max, z_min, z_range,  skewness, kurtosis, n #,sigma 
      del slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4

      # remove rows with any NaNs
      towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

      # make a header string for the output file
      header = 'x, y, z_mean, z_max, z_min, z_range, sigma, skewness, kurtosis, n, slope, intercept, r_value, p_value, std_err, d, l, wmax, wmean, rms1, rms2, Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4'


   #==============================================================================
   print "(4) Writing data to file ..."

   # create a string for the output file
   outfile = infile+'_zstat_detrend'+str(detrend)+'_outres'+str(out)+'_proctype'+str(proctype)+'_mxpts'+str(mxpts)+'_minpts'+str(minpts)+'.xyz' 

   # write the data to the file
   pysesa.write.txtwrite(outfile, towrite, header)

   # stop the clock
   if os.name=='posix': # true if linux/mac
      elapsed = (time() - start1)
   else: # windows
      elapsed = (clock() - start1)

   print "Done! %s points decimated to %s points. Program ran for %s seconds" % (str(orig_pts), str(len(towrite)), str(elapsed))

# =========================================================
# ========================== parsing ======================
# =========================================================

## get list of input arguments and pre-allocate arrays
#argv = sys.argv[1:]
#infile = ''; out = ''; detrend = ''; proctype=''
#mxpts = ''; res = ''; nbin = ''; lentype = ''; 
#minpts = ''; taper = ''; prc_overlap = ''

## parse inputs to variables
#try:
#   opts, args = getopt.getopt(argv,"hi:o:d:x:M:r:b:l:t:v:")
#except getopt.GetoptError:
#      print 'py_sesa.py -i <infile> -o <output size> -d <type of detrend> -x <processing type> -M <max. pts> -m <min. pts> -r <res for gridding> -b <number of bins> -l <lengthscale type 1=l<0.5, 0=l<0> -t <taper type 1=hanning 2=hamming 3=blackman 4=bartlett> -v <percent overlap>'
#      sys.exit(2)
#for opt, arg in opts:
#   if opt == '-h':
#      print 'py_sesa.py -i <infile> -o <output size> -d <type of detrend> -x <processing type> -M <max. pts> -m <min. pts> -r <res for gridding> -b <number of bins> -l <lengthscale type 1=l<0.5, 0=l<0> -t <taper type 1=hanning 2=hamming 3=blackman 4=bartlett> -v <percent overlap>'
#      sys.exit()
#   elif opt in ("-i"):
#      infile = arg
#   elif opt in ("-d"):
#      detrend = arg
#   elif opt in ("-o"):
#      out = arg
#   elif opt in ("-x"):
#      proctype = arg
#   elif opt in ("-M"):
#      mxpts = arg
#   elif opt in ("-m"):
#      minpts = arg
#   elif opt in ("-r"):
#      res = arg
#   elif opt in ("-b"):
#      nbin = arg
#   elif opt in ("-l"):
#      lentype = arg
#   elif opt in ("-t"):
#      taper = arg
#   elif opt in ("-v"):
#      prc_overlap = arg
                        

