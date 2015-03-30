"""
                      __            
 _      ______  _____/ /_____  _____
| | /| / / __ \/ ___/ //_/ _ \/ ___/
| |/ |/ / /_/ / /  / ,< /  __/ /    
|__/|__/\____/_/  /_/|_|\___/_/     
                                    
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

from scipy.interpolate import griddata

import RunningStats
import spec

import warnings
warnings.filterwarnings("ignore")

# =========================================================
cdef class worker:
   """
   Returns an instance. All spatial and spectral parameters requested through .getdata()
   """
   cdef object data, est, myOdr, myModel

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float32_t, ndim=2] points, int proctype, float out, int order=1, float res=0.25, str method='nearest', int nbin=20, int lentype=1, int taper=1): 

      # pre-allocate arrays
      cdef list filled
      cdef list ft

      points = points[np.where(np.logical_not(np.any(np.isnan(points),axis=1)))[0],:]

      # get the (x,y) centroid
      cdef float x = np.min(points[:,0]) + ( np.max(points[:,0]) - np.min(points[:,0]) )/2      
      cdef float y = np.min(points[:,1]) + ( np.max(points[:,1]) - np.min(points[:,1]) )/2                   
                    
      # detrend the x and y to help out the floating point arithmetic    
      points[:,:2] = points[:,:2] - np.min(points[:,:2],axis=0)

      # pre-allocate some more arrays and get the ranges of x and y
      cdef float nx = np.max(points[:,0]) - np.min(points[:,0])
      cdef float ny = np.max(points[:,1]) - np.min(points[:,1])

      cdef float lenx

      # enfore square matrix in gridding
      if ny!=nx:
         lenx = np.ceil(nx/(res*(nx/ny)))
      else:
         lenx = np.ceil(nx/res)
         
      if np.isnan(lenx):
         if proctype<3:
            self.data = np.zeros(35).tolist()
            return
         else:
            self.data = np.zeros(15).tolist()
            return

      # now we can pre-allocate the rest of the arrays
      cdef np.ndarray[np.float32_t, ndim=2] dat = np.empty((lenx,lenx), dtype=np.float32)
      cdef np.ndarray[np.float64_t, ndim=2] Xr = np.empty((len(points),3), dtype=np.float64)      
      cdef np.ndarray[np.float64_t, ndim=2] grid_x = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] grid_y = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] Zf = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] dt = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] resid = np.empty((len(points),), dtype=np.float64)

      # do the gridding
      if ny>nx:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res*(nx/ny)), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )      
      elif ny<nx:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res*(ny/nx)) )      
      else:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )  

      dat = griddata(points[:,:2], points[:,2], (grid_x, grid_y), method=method) #.astype('float32')

      # order=0, no plane detrending is carried out 
      if order==0:

         dt = (dat - np.nanmean(dat.flatten())).astype('float64')
         
         rs1 = RunningStats.RunningStats()
         # global stats, not detrended
         for k in points[:,2].astype('float64'):
            rs1.Push(k)
 
         rs2 = RunningStats.RunningStats()
         # global stats, detrended
         for k in dt.flatten(): #.astype('float64'):
            rs2.Push(k)

      # order=1, OLS plane detrending is carried out                               
      elif order==1:

         import statsmodels.api as sm
      
         ## fit a OLS model with intercept on x and y
         Xr = sm.add_constant(points[:,:2])
         est = sm.OLS(points[:,2], Xr).fit()    

         # plot the hyperplane by evaluating the parameters on the grid
         resid = est.resid
         Zf = self._func(est.params, [grid_x, grid_y])

         dt = (dat-Zf)#.astype('float32')
         
         rs1 = RunningStats.RunningStats()
         # global stats, not detrended
         for k in points[:,2].astype('float64'):
            rs1.Push(k)
 
         rs2 = RunningStats.RunningStats()
         # global stats, detrended
         for k in resid: #.astype('float64'):
            rs2.Push(k)
                 
      # order=2, RLM plane detrending is carried out                                       
      elif order==2:

         import statsmodels.api as sm
      
         ## fit a OLS model with intercept on x and y
         Xr = sm.add_constant(points[:,:2])
         est = sm.RLM(points[:,2], Xr).fit()    

         # plot the hyperplane by evaluating the parameters on the grid
         Zf = self._func(est.params, [grid_x, grid_y])
         resid = est.resid

         dt = (dat-Zf)#.astype('float32')
         
         rs1 = RunningStats.RunningStats()
         # global stats, not detrended
         for k in points[:,2].astype('float64'):
            rs1.Push(k)
 
         rs2 = RunningStats.RunningStats()
         # global stats, detrended
         for k in resid: #.astype('float64'):
            rs2.Push(k)
         
      # order=3, ODR plane detrending is carried out                                       
      elif order==3:

         from scipy.odr import Model, ODR, Data

         myModel = Model(self._func)
         myOdr = ODR(Data([points[:,0], points[:,1]], points[:,2]), myModel, beta0 = [1,1,1])
         myOdr.set_job(fit_type=2)
         est = myOdr.run()

         myModel = Model(self._func)
         myOdr = ODR(Data([points[:,0], points[:,1]], points[:,2]), myModel, beta0 = est.beta)
         myOdr.set_job(fit_type=0)
         est = myOdr.run()

         # plot the hyperplane by evaluating the parameters on the grid
         Zf = self._func(est.beta, [grid_x, grid_y])

         resid = est.delta[0] + est.delta[1] #np.sqrt(est.delta[0]**2 + est.delta[1]**2)

         dt = (dat-Zf) #.astype('float32')

         rs1 = RunningStats.RunningStats()
         # global stats, not detrended
         for k in points[:,2].astype('float64'):
            rs1.Push(k)
 
         rs2 = RunningStats.RunningStats()
         # global stats, detrended
         for k in resid: #.astype('float64'): #flatten():
            rs2.Push(k)
       
      # otherwise, Savitsky-Golay order 0 detrending is carried out                                       
      else:

         import sgolay
               
         Zf = sgolay.sgolay2d( dat, 3, order=0).getdata() 

         Zf[Zf>(np.mean(Zf)+3*np.std(Zf))] = np.nan
         Zf[Zf<(np.mean(Zf)-3*np.std(Zf))] = np.nan

         Zf[np.isnan(Zf)] = np.nanmean(Zf)

         dt = (dat-Zf) #.astype('float32')
         resid = np.random.choice(dt.flatten(),len(points)).astype('float64')

         rs1 = RunningStats.RunningStats()
         # global stats, not detrended
         for k in points[:,2].astype('float64'):
            rs1.Push(k)
 
         rs2 = RunningStats.RunningStats()
         # global stats, detrended
         for k in resid.flatten():#.astype('float64'):
            rs2.Push(k)

      # calls the spectral analysis routine
      ft = spec.spec(dt, nbin, res, proctype, lentype, taper).getdata()    #.astype('float64')  
    
      # compile all parameters into a list
      filled = [x, y, rs1.Mean(), np.max(points[:,2]), np.min(points[:,2]), np.max(points[:,2])-np.min(points[:,2]), rs1.StandardDeviation(), rs2.StandardDeviation(), rs1.Skewness(), rs2.Skewness(), rs1.Kurtosis(), rs2.Kurtosis(), rs2.Mean(), len(points)] + ft

      self.data = filled

      rs1.Clear()
      rs2.Clear()

      return

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   def _func(self,beta,data):
      '''
      minimization function for plane fitting
      '''
      x,y = data
      a,b,c = beta
      return a*x+b*y+c

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      return self.data


