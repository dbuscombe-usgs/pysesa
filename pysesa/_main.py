"""
    ____  __  __________  _________ _
   / __ \/ / / / ___/ _ \/ ___/ __ `/
  / /_/ / /_/ (__  )  __(__  ) /_/ / 
 / .___/\__, /____/\___/____/\__,_/  
/_/    /____/   

+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
  _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _  
 / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ 
( d | b | u | s | c | o | m | b | e | @ | u | s | g | s | . | g | o | v )
 \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ 

+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

"""

# =========================================================
# =============== import libraries ======================
# =========================================================

from __future__ import division
import numpy as np
from joblib import Parallel, delayed, cpu_count
from time import clock, time
import os, sys #, getopt

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

import worker
import partition

import warnings
warnings.filterwarnings("ignore")

# =========================================================
# =============== begin subfunctions ======================
# =========================================================

# =========================================================
def ascol( arr ):
   '''
   reshapes row matrix to be a column matrix (N,1).
   '''
   if len( arr.shape ) == 1: arr = arr.reshape( ( arr.shape[0], 1 ) )
   return arr

#==================================================
def get_stats(pts, proctype, out, order, res, method, nbin, lentype, taper):
   '''
   call the analysis routine. Gets called by the parallel processing queue
   '''
   return worker.worker(pts, proctype, out, order, res, method, nbin, lentype, taper).getdata()

#==================================================        
def txtread(infile):
   '''
   custom fast(er than numpy's genfromtxt) txt file to numpy array
   '''
   f = open(infile, 'rb'); data = f.read(); f.close()
   data = data.splitlines()
   if len(data[0].split(','))>1: # the file is comma delimited
      return np.array([x.split(',',2)[0:3] for x in data], dtype='float32')
   else:
      try: # space delimited
         return np.array([x.split(' ',2)[0:3] for x in data], dtype='float32')
      except: #tab delimited
         return np.array([x.split('\t',2)[0:3] for x in data], dtype='float32')


# =========================================================
# =============== start inputs ======================
# =========================================================
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
    ..   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _  
    ..  / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ 
    .. ( d | b | u | s | c | o | m | b | e | @ | u | s | g | s | . | g | o | v )
    ..  \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ 
    .. 
    .. +-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
    .. |U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
    .. +-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
    .. 
"""

__all__ = [
    'pysesa',
    'ascol',
    'get_stats',
    'txtread',
    ]                
                   
#################################################
def pysesa(infile, out, order, proctype, mxpts, res, nbin, lentype, taper, prc_overlap):  
             
   # if no input file given, exit
   if not infile:
      print 'An input file is required'
      sys.exit(2)

   print 'Input file is %s' % (infile)

   if prc_overlap:
      prc_overlap = np.asarray(prc_overlap,float)
      print 'Percent overlap is %s' % (str(prc_overlap))

   if out:
      out = np.asarray(out,float)
      print 'Output grid size is %s' % (str(out))

   if order:
      order = np.asarray(order,int)
      print 'Order for detrend is %s' % (str(order))

   if proctype:
      proctype = np.asarray(proctype,int)
      if proctype==0:
         proctype = 1

      if proctype==1:
         print 'focal stats, lengthscale and spectral parameters (no smoothing)'
      elif proctype==2:
         print 'focal stats, lengthscale and spectral parameters (with smoothing)'
      elif proctype==3:
         print 'focal stats and lengthscale'

   if res:
      res = np.asarray(res,float)
      print 'Res. is %s' % (str(res))

   if mxpts:
      mxpts = np.asarray(mxpts,int)
      print 'Max points per window is %s' % (str(mxpts))

   if nbin:
      nbin = np.asarray(nbin,float)
      print 'Number of bins is %s' % (str(nbin))

   if lentype:
      lentype = int(lentype)
      if lentype==0:
         print "lengthscale type: l<0"
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
      out = 1
      print '[Default] Output grid size is %s' % (str(out))

   if not order and order!=0:
      order = 3
      print '[Default] Order for detrend is %s (ODR plane)' % (str(order))

   if not proctype:
      proctype = 1
      print '[Default] Processing focal stats, lengthscale and spectral parameters (no smoothing)'

   if not res:
      res = 0.05
      print '[Default] Res. is %s' % (str(res))

   if not mxpts:
      mxpts = 256
      print '[Default] Max points per window is %s' % (str(mxpts))

   if not nbin:
      nbin = 20
      print '[Default] Number of bins is %s' % (str(nbin))

   if not lentype:
      if lentype != 0:
         lentype = 1
         print "[Default] lengthscale type: l>0.5"
 
   if not taper:
      taper = 1
      print '[Default] Hanning taper'


   # max distance for nearest neighbour search
   # specified by spacing (output res) and prc overlap
   win =  np.multiply(out/2, 1+(prc_overlap/100))

   method = 'nearest'

   # start timer
   if os.name=='posix': # true if linux/mac or cygwin on windows
      start1 = time()
   else: # windows
      start1 = clock()

   #==============================================================================
   print "(1) Reading data from file ..."
   #toproc = np.genfromtxt(infile)
   #toproc = txtread(infile)
   toproc = pysesa_read.txtread(infile)

   orig_pts = len(toproc)
   
   nr_pts = pysesa_partition.partition(toproc, out, res, mxpts, win).getdata()

   # start 2nd timer
   if os.name=='posix': # true if linux/mac or cygwin on windows
      start2 = time()
   else: # windows
      start2 = clock()

   #==============================================================================
   print "(3) Processing in parallel using %s processors ... " % (str(cpu_count()))

   try: #parallel processing with all available cores
      w = Parallel(n_jobs=cpu_count(), verbose=0)(delayed(get_stats)(toproc[nr_pts[k],:3], proctype, out, order, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts))) #backend="threading"
   except: #fall back to serial
      w = Parallel(n_jobs=1, verbose=0)(delayed(get_stats)(toproc[nr_pts[k],:3], proctype, out, order, res, method, nbin, lentype, taper) for k in xrange(len(nr_pts)))

   # parse to variables
   if proctype==1: 
      x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_d,ku,ku_d,zmean_d,n,w2,gamma,r_value,p_value,std_err,l,rms1,rms2,wav,Z,E,sigma,T0_1,T0_2,sw1,sw2,m0,m1,m2,m3,m4 = zip(*w)

      w2 = np.asarray(w2)
      gamma = np.asarray(gamma)
      w2 = (10**w2)**(1/(2-gamma))
      
   elif proctype==2: 
      x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_d,ku,ku_d,zmean_d,n,w2,gamma,r_value,p_value,std_err,l,rms1,rms2,wav,Z,E,sigma,T0_1,T0_2,sw1,sw2,m0,m1,m2,m3,m4 = zip(*w)

      w2 = np.asarray(w2)
      gamma = np.asarray(gamma)
      w2 = (10**w2)**(1/(2-gamma))
      
   elif proctype==3: 
      x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_d,ku,ku_d,zmean_d,n,l = zip(*w)  

   del w
   toproc = None
   nr_pts = None

   # combine into single matrix for writing to file
   if (proctype==1) or (proctype==2):
      towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)),ascol(np.asarray(zmean)),ascol(np.asarray(zmax)),ascol(np.asarray(zmin)),ascol(np.asarray(range)),ascol(np.asarray(stdev)),ascol(np.asarray(stdev_d)),ascol(np.asarray(sk)),ascol(np.asarray(sk_d)),ascol(np.asarray(ku)),ascol(np.asarray(ku_d)),ascol(np.asarray(zmean_d)),ascol(np.asarray(w2)),ascol(np.asarray(gamma)),ascol(np.asarray(r_value)),ascol(np.asarray(p_value)),ascol(np.asarray(std_err)),ascol(np.asarray(l)),ascol(np.asarray(rms1)),ascol(np.asarray(rms2)),ascol(np.asarray(wav)),ascol(np.asarray(Z)),ascol(np.asarray(E)),ascol(np.asarray(sigma)),ascol(np.asarray(T0_1)),ascol(np.asarray(T0_2)),ascol(np.asarray(sw1)),ascol(np.asarray(sw2)),ascol(np.asarray(m0)),ascol(np.asarray(m1)),ascol(np.asarray(m2)),ascol(np.asarray(m3)),ascol(np.asarray(m4)),ascol(np.asarray(n)) ))
      del x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_d,ku,ku_d,zmean_d,w2,gamma,r_value,p_value,std_err,l,rms1,rms2,wav,n,Z,E,sigma,T0_1,T0_2,sw1,sw2,m0,m1,m2,m3,m4 

   elif proctype==3:
      towrite = np.hstack(( ascol(np.asarray(x)),ascol(np.asarray(y)),ascol(np.asarray(zmean)),ascol(np.asarray(zmax)),ascol(np.asarray(zmin)),ascol(np.asarray(range)),ascol(np.asarray(stdev)),ascol(np.asarray(stdev_d)),ascol(np.asarray(sk)),ascol(np.asarray(sk_d)),ascol(np.asarray(ku)),ascol(np.asarray(ku_d)),ascol(np.asarray(zmean_d)),ascol(np.asarray(l)), ascol(np.asarray(n)) ))
      del x,y,zmean,zmax,zmin,range,stdev,stdev_d,sk,sk_d,ku,ku_d,zmean_d,l,n

   # remove rows with any NaNs
   towrite = towrite[np.where(np.logical_not(np.any(np.isnan(towrite),axis=1)))[0],:]

   # stop the clock
   if os.name=='posix': # true if linux/mac
      elapsed2 = (time() - start2)
   else: # windows
      elapsed2 = (clock() - start2)

   print "... processing took %s seconds" % (str(elapsed2))

   # write all of this to file
   #==============================================================================
   outfile = infile+'_zstat_order'+str(order)+'_outres'+str(out)+'_proctype'+str(proctype)+'_lentype'+str(lentype)+'_mxpts'+str(mxpts)+'_nbin'+str(nbin)+'.xyz' 
   print "(4) Writing processed points to file ... %s" % (outfile)
   # write points out to a file formated
   if (proctype==1) or (proctype==2):
      with open(outfile, 'w') as f:
         np.savetxt(f, towrite[np.where(towrite[:,-1])[0],:], delimiter=' ', fmt="%8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f")
   elif proctype==3:
      with open(outfile, 'w') as f:
         np.savetxt(f, towrite[np.where(towrite[:,-1])[0],:], delimiter=' ', fmt="%8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f")

   # stop the clock
   if os.name=='posix': # true if linux/mac
      elapsed = (time() - start1)
   else: # windows
      elapsed = (clock() - start1)

   print "Done! %s points decimated to %s points. Program ran for %s seconds" % (str(orig_pts), str(len(towrite)), str(elapsed))


#   toproc = toproc[~np.isnan(toproc).any(axis=1)]

#   xmin = np.min(toproc[:,0])
#   xmax = np.max(toproc[:,0])
#   ymin = np.min(toproc[:,1])
#   ymax = np.max(toproc[:,1])

#   orig_pts = len(toproc)

#   #==============================================================================
#   print "(2) KD-tree for nearest-neighbour lookup ..."
#   # make vector of points decimated to 'out' metres
#   # to be used as search nodes
#   x = np.arange(xmin, xmax, out)
#   y = np.arange(ymin, ymax, out)
#   xx, yy = np.meshgrid(x, y)
#   p = list(np.vstack([xx.flatten(),yy.flatten()]).transpose())

#   # format points for kd-tree
#   allpoints = zip(toproc[:,0].ravel(), toproc[:,1].ravel())

#   # find all points within 'out' metres of each centroid in p 
#   xvec = np.arange(xmin-2*res,xmax+2*res)
#   yvec = np.arange(ymin-2*res,ymax+2*res)
#   nr_pts = partition.partition(toproc, allpoints, p, xvec, yvec, out, res, mxpts, win).getdata()

#   del allpoints, x, y, xx, yy, p

# for debug
#order = 3
#res = 0.05
#proctype = 1
#out = 0.5
#infile='example_data.txt'
#mxpts = 256
#nbin = 20
#lentype = 1 # 1 = l<0.5, else l<0
#taper = 1
#prc_overlap = 0

## get list of input arguments and pre-allocate arrays
#argv = sys.argv[1:]
#infile = ''; out = ''; order = ''; proctype=''
#mxpts = ''; res = ''; nbin = ''; lentype = ''; 
#taper = ''; prc_overlap = ''

## parse inputs to variables
#try:
#   opts, args = getopt.getopt(argv,"hi:o:p:x:m:r:b:l:t:v:")
#except getopt.GetoptError:
#      print 'py_sesa.py -i <direc> -o <output size> -p <order of detrend> -x <processing type> -m <max. pts> -r <res for gridding> -b <number of bins> -l <lengthscale type 1=l<0.5, 0=l<0> -t <taper type 1=hanning 2=hamming 3=blackman 4=bartlett> -v <percent overlap>'
#      sys.exit(2)
#for opt, arg in opts:
#   if opt == '-h':
#      print 'py_sesa.py -i <direc> -o <output size> -p <order of detrend> -x <processing type> -m <max. pts> -r <res for gridding> -b <number of bins> -l <lengthscale type 1=l<0.5, 0=l<0> -t <taper type 1=hanning 2=hamming 3=blackman 4=bartlett> -v <percent overlap>'
#      sys.exit()
#   elif opt in ("-i"):
#      infile = arg
#   elif opt in ("-p"):
#      order = arg
#   elif opt in ("-o"):
#      out = arg
#   elif opt in ("-x"):
#      proctype = arg
#   elif opt in ("-m"):
#      mxpts = arg
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
