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

#"""                       
# ___      ___ ___ ___   _     _   _ 
#| _ \_  _/ __| __/ __| /_\   (_) (_)
#|  _/ || \__ \ _|\__ \/ _ \   _   _ 
#|_|  \_, |___/___|___/_/ \_\ (_) (_)
#     |__/                           
#                     __  _       __
#   _________  ____ _/ /_(_)___ _/ /
#  / ___/ __ \/ __ `/ __/ / __ `/ / 
# (__  ) /_/ / /_/ / /_/ / /_/ / /  
#/____/ .___/\__,_/\__/_/\__,_/_/   
#    /_/                            

#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
#|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
#"""

# import libraries
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
import RunningStats

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
      rs1 = RunningStats.RunningStats()
      for k in points[:,2]:
         rs1.Push(k)
      #z_mean,z_max,z_min,z_range,sigma,skewness,kurtosis,n
      return [rs1.Mean(), np.max(points[:,2]), np.min(points[:,2]), np.max(points[:,2])-np.min(points[:,2]), rs1.StandardDeviation(), rs1.Skewness(), rs1.Kurtosis(), len(points)] 



