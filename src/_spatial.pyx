# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8


# import libraries
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

from scipy.stats import kurtosis
from scipy.stats import skew

# suppress divide and invalid warnings
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

# =========================================================
cdef class spatial:
   '''
   Calculate spatial statistics of a Nx3 point cloud

   Syntax
   ----------
   stats = pysesa.spatial(points).getdata()

   centroids = pysesa.spatial(points).getcentroid()

   stats = pysesa.spatial(points).getstats()

   Parameters
   -------------
   points : ndarray
   	Nx3 point cloud

   Returns [requested through .getdata()]
   ---------------------------------------
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
   --------------------------------------------
   self.centroid: list
   	1x3 point cloud centroid [x,y,z]

   Returns [requested through .getstats()]
   ----------------------------------------
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

   cdef object data, rs1, centroid, stats

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] points):
      '''
      Calculate spatial statistics of a Nx3 point cloud

      Syntax
      ----------
      stats = pysesa.spatial(points).getdata()

      centroids = pysesa.spatial(points).getcentroid()

      stats = pysesa.spatial(points).getstats()

      Parameters
      -------------
      points : ndarray
   	   Nx3 point cloud

      Returns [requested through .getdata()]
      -----------------------------------------
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
      --------------------------------------------
      self.centroid: list
   	   1x3 point cloud centroid [x,y,z]

      Returns [requested through .getstats()]
      -----------------------------------------
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

      cdef float x, y, z_mean,z_maz,z_min,z_range,sigma,skewness,kurtosis,n

      x,y,z = self._get_centroid(points)

      z_mean,z_max,z_min,z_range,sigma,skewness,kurtosis,n = self._get_stats(points)

      self.stats = [z_mean,z_max,z_min,z_range,sigma,skewness,kurtosis,n]
      self.centroid = [x,y,z]

      #x,y,z_mean,z_max,z_min,z_range,sigma,skewness,kurtosis,n
      self.data = [x,y,z_mean,z_max,z_min,z_range,sigma,skewness,kurtosis,n]

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      '''
      Calculate spatial statistics of a Nx3 point cloud

      Syntax
      ----------
      stats = pysesa.spatial.getdata()

      Parameters
      -------------
      self : instance
   	   pysesa.spatial instance

      Returns [requested through .getdata()]
      ---------------------------------------
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

      '''
      return self.data

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getcentroid(self):
      '''
      Calculate spatial statistics of a Nx3 point cloud

      Syntax
      ----------
      centroids = pysesa.spatial.getcentroid()

      Parameters
      -------------
      self : instance
   	   pysesa.spatial instance

      Returns [requested through .getcentroid()]
      --------------------------------------------
      self.centroid: list
   	   1x3 point cloud centroid [x,y,z]

      '''
      return self.centroid

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getstats(self):
      '''
      Calculate spatial statistics of a Nx3 point cloud

      Syntax
      ----------
      stats = pysesa.spatial.getstats()

      Parameters
      ------------
      self : instance
   	   pysesa.spatial instance

      Returns [requested through .getstats()]
      -----------------------------------------
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
      return self.stats

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _get_centroid(self, np.ndarray[np.float64_t, ndim=2] points):
      '''
      Calculate centroid of a Nx3 point cloud

      Syntax
      ----------
      centroids = pysesa.spatial.getcentroid(points)

      Parameters
      ------------
      points : ndarray
   	   Nx3 point cloud

      Returns [requested through .getcentroid()]
      --------------------------------------------
      self.centroid: list
   	   1x3 point cloud centroid [x,y,z]

      '''
      cdef float x = np.min(points[:,0]) + ( np.max(points[:,0]) - np.min(points[:,0]) )/2
      cdef float y = np.min(points[:,1]) + ( np.max(points[:,1]) - np.min(points[:,1]) )/2
      cdef float z = np.min(points[:,2]) + ( np.max(points[:,2]) - np.min(points[:,2]) )/2

      return [x, y, z]

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _get_stats(self, np.ndarray[np.float64_t, ndim=2] points):
      '''
      Calculate spatial statistics of a Nx3 point cloud

      Syntax
      ----------
      stats = pysesa.spatial.getstats(points)

      Parameters
      ------------
      points : ndarray
   	   Nx3 point cloud

      Returns [requested through .getstats()]
      ----------------------------------------
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

      return [np.mean(points[:,2]), np.max(points[:,2]), np.min(points[:,2]), np.max(points[:,2])-np.min(points[:,2]), np.std(points[:,2]), skew(points[:,2]), kurtosis(points[:,2]), len(points)]
