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

##"""
## ___      ___ ___ ___   _     _   _ 
##| _ \_  _/ __| __/ __| /_\   (_) (_)
##|  _/ || \__ \ _|\__ \/ _ \   _   _ 
##|_|  \_, |___/___|___/_/ \_\ (_) (_)
##     |__/                           
##                      __  _ __  _           
##    ____  ____ ______/ /_(_) /_(_)___  ____ 
##   / __ \/ __ `/ ___/ __/ / __/ / __ \/ __ \
##  / /_/ / /_/ / /  / /_/ / /_/ / /_/ / / / /
## / .___/\__,_/_/   \__/_/\__/_/\____/_/ /_/ 
##/_/                                         

##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
##|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
##+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
##|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

##"""

# import libraries
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

try:
   from pykdtree.kdtree import KDTree
   pykdtree=1   
except:
   #print "install pykdtree for faster kd-tree operations: https://github.com/storpipfugl/pykdtree"
   from scipy.spatial import cKDTree as KDTree
   pykdtree=0   

# =========================================================
cdef class partition:

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
   -----------
   toproc : ndarray
   	Nx3 point cloud

   Other Parameters
   -----------------
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
   
   cdef object data, mytree, tree, tree2

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   #==================================================
   def __init__(self, np.ndarray[np.float32_t, ndim=2] toproc, float out=0.5, float mxpts=1024, float minpts=16, float prc_overlap=0): #float res=0.05, int bp=0):

      '''
      Partition a Nx3 point cloud into M windows of nx3 points
      with specified spacing between centroids of adjacent windows
      and with specified overlap between windows.
      Implemented using a binary search tree for fast nearest neighbour 
      point check with optional boundary pruning

      Syntax
      ----------
      nr_pts = pysesa.partition(toproc, out, mxpts, minpts, prc_overlap).getdata()

      Parameters
      ----------
      toproc : ndarray
   	   Nx3 point cloud

      Other Parameters
      ----------
      out : float, *optional* [default = 0.5]
   	   output grid resolution
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

      #get size of window
      cdef float win =  (out/2) * (1+(prc_overlap/100)) #np.multiply(out/2, 1+(prc_overlap/100))

      #remove any rows containing NaNs
      toproc = toproc[~np.isnan(toproc).any(axis=1)]

      # pre-allocate arrays     
      cdef float xmin = np.min(toproc[:,0])
      cdef float xmax = np.max(toproc[:,0])
      cdef float ymin = np.min(toproc[:,1])
      cdef float ymax = np.max(toproc[:,1])
      cdef int orig_pts = len(toproc)
      cdef int lenx = np.ceil((xmax-xmin)/out)
      cdef int leny = np.ceil((ymax-ymin)/out)
      cdef int lenx2 = lenx/2    
      cdef int leny2 = leny/2   
            
      cdef np.ndarray[np.float64_t, ndim=1] xvec = np.empty((lenx2,), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] yvec = np.empty((leny2,), dtype=np.float64)   

#      cdef np.ndarray[np.float32_t, ndim=2] p = np.empty((lenx*leny,2), dtype=np.float32)   

      cdef np.ndarray[np.float64_t, ndim=2] xx = np.empty((leny,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] yy = np.empty((leny,lenx), dtype=np.float64)

      # create output grid
      x = np.arange(xmin, xmax, out)
      y = np.arange(ymin, ymax, out)
      xx, yy = np.meshgrid(x, y)
      p = np.vstack([xx.flatten(),yy.flatten()]).transpose().astype('float32')
    
      # find all points within 'out' metres of each centroid in p 
      xvec = np.arange(xmin-2*out,xmax+2*out)
      yvec = np.arange(ymin-2*out,ymax+2*out)            

      # pre-allocate more arrays
      cdef int k
      #cdef tuple cx, cy
      cdef list indices_list
      #cdef list w = []
      #cdef list indices2 = []

      cdef np.ndarray[np.float64_t, ndim=2] xp = np.empty((len(yvec), len(xvec)), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] yp = np.empty((len(yvec), len(xvec)), dtype=np.float64)
#      cdef np.ndarray[np.float32_t, ndim=1] dist2 = np.empty((len(xvec)*len(yvec),), dtype=np.float32)
# 
      mytree = KDTree(toproc[:,:2]) #, leafsize=len(toproc)/100)
      if pykdtree==1:
         complete=0
         while complete==0:
            try:
               dist, indices = mytree.query(p,mxpts, distance_upper_bound=win)
               if 'indices' in locals():
                  complete=1
            except:
               mxpts = np.max([1,mxpts-2])
               print "memory error, using %s max points" % (str(mxpts))         

      else:
         try:
            dist, indices = mytree.query(p,mxpts, distance_upper_bound=win, n_jobs=-1)
            #del p
         except:
            complete=0
            while complete==0:
               try:
                  dist, indices = mytree.query(p,mxpts, distance_upper_bound=win)
               except:
                  mxpts = np.max([1,mxpts-2])
                  print "memory error, using %s max points" % (str(mxpts))         

         finally:
            import dask.array as da
            complete=0
            while complete==0:
               try:
                  #dask implementation
                  dat = da.from_array(toproc[:,:2].astype('float64'), chunks=1000)
                  mytree = KDTree(dat, leafsize=len(dat)/100) # adding this leafsize option speeds up considerably
                  dbp = da.from_array(np.asarray(p), chunks=1000) 
                  #del p  
                  dist, indices = mytree.query(dbp.compute(),mxpts, distance_upper_bound=win) #.astype('float32')
                  del dbp
                  if 'indices' in locals():
                     complete=1
               except:
                  mxpts = np.max([1,mxpts-2])
                  print "memory error, using %s max points" % (str(mxpts))         


      # remove any indices associated with 'inf' distance
      indices = np.squeeze(indices[np.where(np.all(np.isinf(dist),axis=1) ==  False),:])
      dist = np.squeeze(dist[np.where(np.all(np.isinf(dist),axis=1) ==  False),:])

      # define null indices
      try:
         indices[indices==len(toproc)] = -999
      except:
         indices[indices==len(dat)] = -999  #dask implementation    

      indices_list = indices.tolist()
      del indices

      # remove null records
      #for k in xrange(len(indices_list)):
      for k from 0 <= k < len(indices_list):
         indices_list[k] = [x for x in indices_list[k] if x!=-999]

      #for k in xrange(len(indices_list)):
      for k from 0 <= k < len(indices_list):
         indices_list[k] = [x for x in indices_list[k] if x<len(toproc)]

      #try:
      #   # get the centroids
      #   #for k in xrange(len(indices_list)):
      #   for k from 0 <= k < len(indices_list):
      #      #w.append(np.mean([allpoints[i] for i in indices_list[k]], axis=0))
      #      w.append(np.mean([toproc[i,:2] for i in indices_list[k]], axis=0))

      #except:
      #   # get the centroids
      #   #for k in xrange(len(indices_list)):
      #   for k from 0 <= k < len(indices_list):
      #      w.append(dat[indices_list[k],:2].mean(axis=0).compute()) #dask implementation

      #cx,cy = zip(*w)
      #del w

      #yp, xp = np.meshgrid(yvec, xvec)

      # Use KDTree to answer the question: "which point of set (x,y) is the
      # nearest neighbors of those in (xp, yp)"

      #if pykdtree==1:
      #   dist2, _ = mytree.query(np.c_[xp.ravel(), yp.ravel()].astype('float32'), k=1)
      #else:
      #   try:      
      #      dist2, _ = mytree.query(np.c_[xp.ravel(), yp.ravel()].astype('float32'), k=1, n_jobs=-1)
      #   except:
      #      dist2, _ = mytree.query(np.c_[xp.ravel(), yp.ravel()].astype('float32'), k=1)      

      #if bp==1: #boundary pruning
      #   # Select points sufficiently far away (use hypoteneuse of the triangle made by res and res)
      #   tree2 = KDTree(np.c_[xp.ravel()[(dist2 > np.hypot(out, out))], yp.ravel()[(dist2 > np.hypot(out, out))]])

      #   if pykdtree==1:     
      #      dist3, _ = tree2.query(np.c_[cx,cy], distance_upper_bound=win) #distance_upper_bound=out)
      #   else: 
      #      try:
      #         dist3, _ = tree2.query(np.c_[cx,cy], distance_upper_bound=win, n_jobs=-1) #distance_upper_bound=out)
      #      except:
      #         dist3, _ = tree2.query(np.c_[cx,cy], distance_upper_bound=win)      
      
      #   m2 = np.where(dist3 < np.hypot(out, out))[0]

      #else:
      #   m2 = np.where(dist2 < np.hypot(out, out))[0]

      #m2 = m2[np.where(m2<len(indices_list))[0]]

      ## do the pruning
      ##for k in xrange(len(m2)):
      #for k from 0 <= k < len(m2):
      #   if len(indices_list[m2[k]])>minpts:
      #      indices2.append(indices_list[m2[k]])

      #if len(indices2)<1:
      #   print "no returned windows: either increase 'out' or 'prc_overlap'"

      self.data = indices_list #indices2

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      '''
      Partition a Nx3 point cloud into M windows of nx3 points
      with specified spacing between centroids of adjacent windows
      and with specified overlap between windows.
      Implemented using a binary search tree for fast nearest neighbour 
      point check with boundary pruning

      Syntax
      ----------
      nr_pts = pysesa.partition.getdata()

      Parameters
      ----------
      self : instance
   	   pysesa.partition instance

      Returns
      ----------
      self.data: list
   	   list of M ndarrays, each containing n indices 
           of original point cloud, toproc, to partition space 
           to create M windows
      '''
      return self.data


      #cdef np.ndarray[np.float64_t, ndim=2] indices = np.empty((len(p),mxpts), dtype=np.int64)
       #cdef np.ndarray[np.float64_t, ndim=2] xx = np.empty((int(np.ceil((ymax-ymin)/out)),int(np.ceil((xmax-xmin)/out))), dtype=np.float64)
      #cdef np.ndarray[np.float64_t, ndim=2] yy = np.empty((int(np.ceil((ymax-ymin)/out)),int(np.ceil((xmax-xmin)/out))), dtype=np.float64)
      #dbp = db.from_sequence(p, npartitions = 1000) #dask bag
      #del p #dask

      # format points for kd-tree
      #allpoints = zip(toproc[:,0].ravel(), toproc[:,1].ravel())
      #del toproc
       
      # get the tree for the entire point cloud
      #mytree = cKDTree(allpoints) 
      
      # largest inscribed square has side length = sqrt(2)*radius
      #dist, indices = mytree.query(p,mxpts, distance_upper_bound=win)

      #dist, indices = mytree.query(dbp.compute(),mxpts, distance_upper_bound=win) #dask implementation 1
      
      #cdef int bp=0 #boundary prining flag
      
      # no idea why this doesnt work
      #cdef np.ndarray[np.float64_t, ndim=1] x = np.empty(lenx, dtype=np.float64)
      #cdef np.ndarray[np.float64_t, ndim=1] y = np.empty(leny, dtype=np.float64)
      #try:
      #   tree = KDTree(allpoints, leafsize=len(allpoints)/100)
      #except:
      #   tree = KDTree(dat, leafsize=len(dat)/100) #dask implementation 2 # leafsize=len(dat)
      #   del dat
      #p = list(np.vstack([xx.flatten(),yy.flatten()]).transpose()) 
      #cdef list p
      #cdef list allpoints 
#from scipy.spatial import cKDTree
#import dask.bag as db
   
#data_type = np.float64
#ctypedef np.float64_t data_type_f64t

#ctypedef double (*metric_ptr)(double[:, ::1], np.intp_t, np.intp_t)
            
