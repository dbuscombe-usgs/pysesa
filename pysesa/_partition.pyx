"""
                      __  _ __  _           
    ____  ____ ______/ /_(_) /_(_)___  ____ 
   / __ \/ __ `/ ___/ __/ / __/ / __ \/ __ \
  / /_/ / /_/ / /  / /_/ / /_/ / /_/ / / / /
 / .___/\__,_/_/   \__/_/\__/_/\____/_/ /_/ 
/_/ 

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

# import libraries
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from scipy.spatial import cKDTree

# =========================================================
cdef class partition:

   '''
   binary search tree for fast nearest neighbour point check
   with boundary pruning
   nr_pts = pysesa_partition.partition(toproc, allpoints, p, xvec, yvec, out, res, mxpts).getdata()
   '''
   
   cdef object data, mytree, tree, tree2

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   #==================================================
   def __init__(self, np.ndarray[np.float32_t, ndim=2] toproc, float out=0.5, float res=0.05, float mxpts=256, float win=100):

      toproc = toproc[~np.isnan(toproc).any(axis=1)]
      cdef list p
      cdef list allpoints
      cdef float xmin = np.min(toproc[:,0])
      cdef float xmax = np.max(toproc[:,0])
      cdef float ymin = np.min(toproc[:,1])
      cdef float ymax = np.max(toproc[:,1])
      cdef int orig_pts = len(toproc)
      cdef int lenx = np.ceil((xmax-xmin)/out)
      cdef int leny = np.ceil((ymax-ymin)/out)
      cdef int lenx2 = lenx/2    
      cdef int leny2 = leny/2   
      
#      cdef np.ndarray[np.float64_t, ndim=1] x = np.empty(lenx, dtype=np.float64)
#      cdef np.ndarray[np.float64_t, ndim=1] y = np.empty(leny, dtype=np.float64)

      cdef np.ndarray[np.float64_t, ndim=1] xvec = np.empty((lenx2,), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] yvec = np.empty((leny2,), dtype=np.float64)   
   
      cdef np.ndarray[np.float64_t, ndim=2] xx = np.empty((int(np.ceil((ymax-ymin)/out)),int(np.ceil((xmax-xmin)/out))), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] yy = np.empty((int(np.ceil((ymax-ymin)/out)),int(np.ceil((xmax-xmin)/out))), dtype=np.float64)

      x = np.arange(xmin, xmax, out)
      y = np.arange(ymin, ymax, out)
      xx, yy = np.meshgrid(x, y)
      p = list(np.vstack([xx.flatten(),yy.flatten()]).transpose())

      # format points for kd-tree
      allpoints = zip(toproc[:,0].ravel(), toproc[:,1].ravel())
        
      # find all points within 'out' metres of each centroid in p 
      xvec = np.arange(xmin-2*res,xmax+2*res)
      yvec = np.arange(ymin-2*res,ymax+2*res)            

      cdef int k
      cdef tuple cx, cy
      cdef list indices_list
      cdef list w = []
      cdef list indices2 = []

      cdef np.ndarray[np.float64_t, ndim=2] dist = np.empty((len(p),mxpts), dtype=np.float64)
      cdef np.ndarray[np.int64_t, ndim=2] indices = np.empty((len(p),mxpts), dtype=np.int64)
      cdef np.ndarray[np.float64_t, ndim=2] xp = np.empty((len(yvec), len(xvec)), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] yp = np.empty((len(yvec), len(xvec)), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] dist2 = np.empty((len(xvec)*len(yvec),), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] dist3 = np.empty((len(p),), dtype=np.float64)
 
      # get the tree for the entire point cloud
      mytree = cKDTree(allpoints) 
   
      # largest inscribed square has side length = sqrt(2)*radius
      dist, indices = mytree.query(p,mxpts, distance_upper_bound=win)
      # remove any indices associated with 'inf' distance
      indices = np.squeeze(indices[np.where(np.all(np.isinf(dist),axis=1) ==  False),:])
      dist = np.squeeze(dist[np.where(np.all(np.isinf(dist),axis=1) ==  False),:])

      # define null indices
      indices[indices==len(allpoints)] = -999
      indices_list = indices.tolist()
      del indices

      # remove null records
      #for k in xrange(len(indices_list)):
      for k from 0 <= k < len(indices_list):
         indices_list[k] = [x for x in indices_list[k] if x!=-999]
      
      # get the centroids
      #for k in xrange(len(indices_list)):
      for k from 0 <= k < len(indices_list):
         w.append(np.mean(toproc[indices_list[k],:2],axis=0))

      cx,cy = zip(*w)
      del w

      yp, xp = np.meshgrid(yvec, xvec)

      # Use KDTree to answer the question: "which point of set (x,y) is the
      # nearest neighbors of those in (xp, yp)"
      tree = cKDTree(allpoints)
      dist2, _ = tree.query(np.c_[xp.ravel(), yp.ravel()], k=1)

      # Select points sufficiently far away (use hypoteneuse of the triangle made by res and res)

      tree2 = cKDTree(np.c_[xp.ravel()[(dist2 > np.hypot(res, res))], yp.ravel()[(dist2 > np.hypot(res, res))]])
      dist3, _ = tree.query(np.c_[cx,cy], distance_upper_bound=win) #distance_upper_bound=out)
      m2 = np.where(dist3 < out**2)[0]

      # do the pruning
      #for k in xrange(len(m2)):
      for k from 0 <= k < len(m2):
         indices2.append(indices_list[m2[k]])

      self.data = indices2

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      return self.data

