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
##       __     __                      __
##  ____/ /__  / /_________  ____  ____/ /
## / __  / _ \/ __/ ___/ _ \/ __ \/ __  / 
##/ /_/ /  __/ /_/ /  /  __/ / / / /_/ /  
##\__,_/\___/\__/_/   \___/_/ /_/\__,_/   
##                                  
##      
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

# suppress divide and invalid warnings
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

# =========================================================
cdef class detrend:
   '''
   Detrend a Nx3 point cloud
   by specified method

   Syntax
   ----------
   detrended_pts = pysesa.detrend(points, proctype, res, method).getdata()

   Parameters
   ----------
   points : ndarray
   	Nx3 point cloud
   proctype : int
   	type of detrending.
        1 = remove mean

        2 = remove Ordinary least squares plane

        3 = remove Robust linear model plane

        4 = remove Orthogonal Distance Regression plane

        5 = remove Savitsky-Golay digital filter, order 1


   Other Parameters
   ----------
   res : float, *optional* [default = 0.05]
   	for proctype==4 only
        spatial grid resolution to create a grid
   method : str, *optional* [default = 'nearest']
   	for proctype==4 only
   	gridding type

   Returns
   ----------
   self.data: ndarray
   	Nx3 detrended point cloud

   '''

   cdef object data, est, myOdr, myModel

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float32_t, ndim=2] points, int proctype, float res=0.05, str method='nearest'): 

      '''
      Detrend a Nx3 point cloud
      by specified method

      Syntax
      ----------
      detrended_pts = pysesa.detrend(points, proctype, res, method).getdata()

      Parameters
      ----------
      points : ndarray
   	   Nx3 point cloud
      proctype : int
   	   type of detrending.
           1 = remove mean
           2 = remove Ordinary least squares plane
           3 = remove Robust linear model plane
           4 = remove Orthogonal Distance Regression plane
           5 = remove Savitsky-Golay digital filter, order 1

      Other Parameters
      ----------
      res : float, *optional* [default = 0.05]
   	   for proctype==4 only
           spatial grid resolution to create a grid
      method : str, *optional* [default = 'nearest']
   	   for proctype==4 only
   	   gridding type

      Returns
      ----------
      self.data: ndarray
   	   Nx3 detrended point cloud

      '''
      points = points[np.where(np.logical_not(np.any(np.isnan(points),axis=1)))[0],:]

      # pre-allocate some more arrays and get the ranges of x and y
      cdef float nx = np.max(points[:,0]) - np.min(points[:,0])
      cdef float ny = np.max(points[:,1]) - np.min(points[:,1])

      cdef np.float64_t mz = np.empty((1,), dtype=np.float64)      
      mz = np.mean(points[:,2]).astype('float64')

      cdef float lenx

      # enfore square matrix in gridding
      if ny!=nx:
         lenx = np.ceil(nx/(res*(nx/ny)))
      else:
         lenx = np.ceil(nx/res)
         
      if np.isnan(lenx):
         self.data = points.astype('float64')
         return 

      # now we can pre-allocate the rest of the arrays
      cdef np.ndarray[np.float32_t, ndim=2] dat = np.empty((lenx,lenx), dtype=np.float32)
      cdef np.ndarray[np.float64_t, ndim=2] Xr = np.empty((len(points),3), dtype=np.float64)      
      cdef np.ndarray[np.float64_t, ndim=2] grid_x = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] grid_y = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] Zf = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] dt = np.empty((len(points),), dtype=np.float64)      

      # =0, no plane detrending is carried out 
      if proctype==1:

         dt = (points[:,2] - np.mean(points[:,2])).astype('float64')

      # =1, OLS plane detrending is carried out                               
      elif proctype==2:

         import statsmodels.api as sm
      
         ## fit a OLS model with intercept on x and y
         Xr = sm.add_constant(points[:,:2])
         est = sm.OLS(points[:,2], Xr).fit()    

         # plot the hyperplane by evaluating the parameters on the grid
         dt = (est.resid).astype('float64')   
                 
      # =2, RLM plane detrending is carried out                                       
      elif proctype==3:

         import statsmodels.api as sm
      
         ## fit a RLM model with intercept on x and y
         Xr = sm.add_constant(points[:,:2])
         est = sm.RLM(points[:,2], Xr).fit()    
         dt = (est.resid).astype('float64')   
         
      # =3, ODR plane detrending is carried out                                       
      elif proctype==4:

         from scipy.odr import Model, ODR, Data

         myModel = Model(self._func)
         myOdr = ODR(Data([points[:,0], points[:,1]], points[:,2]), myModel, beta0 = [1,1,1])
         myOdr.set_job(fit_type=2)
         est = myOdr.run()

         myModel = Model(self._func)
         myOdr = ODR(Data([points[:,0], points[:,1]], points[:,2]), myModel, beta0 = est.beta)
         myOdr.set_job(fit_type=0)
         est = myOdr.run()

         dt = -(est.delta[0] + est.delta[1]).astype('float64')   
       
      # otherwise, Savitsky-Golay order 0 detrending is carried out                                       
      else:

         import sgolay       
         from scipy.interpolate import griddata

         # do the gridding
         if ny>nx:
            grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res*(nx/ny)), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )      
         elif ny<nx:
            grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res*(ny/nx)) )      
         else:
            grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )  

         dat = griddata(points[:,:2], points[:,2], (grid_x, grid_y), method=method)
     
         Zf = sgolay.sgolay( dat, 3, order=0).getdata() 

         Zf[Zf>(np.mean(Zf)+3*np.std(Zf))] = np.nan
         Zf[Zf<(np.mean(Zf)-3*np.std(Zf))] = np.nan

         Zf[np.isnan(Zf)] = np.mean(Zf)

         Zf = (dat-Zf)
         dt = np.random.choice(Zf.flatten(),len(points)).astype('float64')   

      points[:,2] = dt + mz

      self.data = points.astype('float64')    

      return

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   def _func(self,beta,data):
      '''
      Minimization function for plane fitting

      Parameters
      ----------
      beta : list
   	   vector of minimization parameters

      data : ndarray
   	   Nx2 coordinates

      Returns
      ----------
      plane: ndarray
   	   Nx2 plane coordinates
      '''
      x,y = data
      a,b,c = beta
      return a*x+b*y+c

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] getdata(self):

      '''
      Detrend a Nx3 point cloud
      by specified method

      Syntax
      ----------
      detrended_pts = pysesa.detrend.getdata()

      Parameters
      ----------
      self : instance
   	   pysesa.detrend instance

      Returns
      ----------
      self.data: ndarray
   	   Nx3 detrended point cloud

      '''
      return self.data



