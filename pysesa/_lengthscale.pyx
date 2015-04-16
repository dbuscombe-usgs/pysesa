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
##    __                 __  __                    __   
##   / /__  ____  ____ _/ /_/ /_  ______________ _/ /__ 
##  / / _ \/ __ \/ __ `/ __/ __ \/ ___/ ___/ __ `/ / _ \
## / /  __/ / / / /_/ / /_/ / / (__  ) /__/ /_/ / /  __/
##/_/\___/_/ /_/\__, /\__/_/ /_/____/\___/\__,_/_/\___/ 
##             /____/                                   

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
from scipy.interpolate import griddata
from scipy.integrate import trapz

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')
import warnings
warnings.filterwarnings("ignore")

# =========================================================
cdef class lengthscale:
   '''
   Calculates the integral lengthscale of a Nx3 point cloud
   using 1 of 3 available methods
   and also returns the tapered 2D grid of 3D pointcloud for spectral analysis

   Syntax
   ----------
   im = pysesa.lengthscale(points, res, lentype, taper, method).getdata()
   lengthscale = pysesa.lengthscale(points, res, lentype, taper, method).getlengthscale()

   Parameters
   ------------
   points : ndarray
   	Nx3 point cloud

   Other Parameters
   ------------------
   res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
   lentype : int, *optional* [default = 0, l<0.5]
   	lengthscale type:
        1, l<0.5

        2, l<1/e

        3, l<0

   taper : int, *optional* [default = Hanning]
   	flag for taper type:
        1, Hanning (Hann)

        2, Hamming

        3, Blackman

        4, Bartlett

   method : str, *optional* [default = 'nearest']
   	gridding type

   Returns [requested through .getdata()]
   -----------------------------------------
   self.data: ndarray
   	tapered 2D grid of 3D pointcloud

   Returns [requested through .getlengthscale()]
   ----------------------------------------------
   self.lengthscale: float
   	integral lengthscale

   '''

   cdef object data, lengthscale

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] points, float res=0.05, int lentype=1, int taper=1, str method='nearest'): 
      '''
      Calculates the integral lengthscale of a Nx3 point cloud
      using 1 of 3 available methods
      and also returns the tapered 2D grid of 3D pointcloud for spectral analysis

      Syntax
      ----------
      im = pysesa.lengthscale(points, res, lentype, taper, method).getdata()
      lengthscale = pysesa.lengthscale(points, res, lentype, taper, method).getlengthscale()

      Parameters
      ----------
      points : ndarray
   	   Nx3 point cloud

      Other Parameters
      ----------
      res : float, *optional* [default = 0.05]
           spatial grid resolution to create a grid
      lentype : int, *optional* [default = 1, l<0.5]
   	   lengthscale type:
           1, l<0.5
           2, l<1/e
           3, l<0
      taper : int, *optional* [default = Hanning]
   	   flag for taper type:
           1, Hanning (Hann)
           2, Hamming
           3, Blackman
           4, Bartlett
      method : str, *optional* [default = 'nearest']
   	   gridding type

      Returns [requested through .getdata()]
      ----------
      self.data: ndarray
   	   tapered 2D grid of 3D pointcloud

      Returns [requested through .getlengthscale()]
      ----------
      self.lengthscale: float
   	   integral lengthscale

      '''

      # pre-allocate some more arrays and get the ranges of x and y
      cdef float nnx = np.max(points[:,0]) - np.min(points[:,0])
      cdef float nny = np.max(points[:,1]) - np.min(points[:,1])
      cdef float lenx
      cdef int nx, ny, i

      # enfore square matrix in gridding
      if nny!=nnx:
         lenx = np.ceil(nnx/(res*(nnx/nny)))
      else:
         lenx = np.ceil(nnx/res)

      cdef np.ndarray[np.float64_t, ndim=2] grid_x = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] grid_y = np.empty((lenx,lenx), dtype=np.float64)

      # do the gridding
      if nny>nnx:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res*(nnx/nny)), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )      
      elif nny<nnx:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res*(nny/nnx)) )      
      else:
         grid_x, grid_y = np.meshgrid( np.arange(np.min(points,axis=0)[0], np.max(points,axis=0)[0], res), np.arange(np.min(points,axis=0)[1], np.max(points,axis=0)[1], res) )  


      im = griddata(points[:,:2], points[:,2], (grid_x, grid_y), method=method)

      im = im - np.nanmean(im)
      ny, nx= np.shape(im)

      #taper
      im = self._taper_im(im, taper)
      self.data = im

      #get integral lengthscale
      out = self._get_lengthscale(im, lentype, res, nx, ny)
      self.lengthscale = out


   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray _radial_data(self, np.ndarray[np.float64_t, ndim=2] data, int annulus_width=1):
       '''
       Efficient radial average of matrix

       Syntax
       ----------
       h = pysesa.lengthscale._radial_data(data, annulus_width)

       Parameters
       ----------
       data : ndarray
   	   2D grid
       annulus_width : int [default = 1]
           angular resolution of sector mean

       Returns
       ----------
       radialdatamean: ndarray
   	   radial average 1D array

       '''
       cdef int npix, npiy, nrad, irad
       cdef double rmax

       npix, npiy = np.shape(data)

       cdef np.ndarray[np.float64_t, ndim=2] r = np.empty((npix, npiy),dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=2] x = np.empty((npix, npiy),dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=2] y = np.empty((npix, npiy),dtype=np.float64)  
       cdef np.ndarray[np.float64_t, ndim=1] x1 = np.empty(npix,dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=1] y1 = np.empty(npiy,dtype=np.float64)

       cdef np.ndarray[np.float64_t, ndim=1] minrad = np.empty(1,dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=1] maxrad = np.empty(1,dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=1] dr = np.empty(1,dtype=np.float64)

       x1 = np.arange(-np.float64(npix/2),np.float64(npix/2))
       y1 = np.arange(-np.float64(npiy/2),np.float64(npiy/2))
       x,y = np.meshgrid(y1,x1)
       r = np.abs(x+1j*y)

       rmax = np.max(r)
       dr = np.abs([x[0,0] - x[0,1]]) * annulus_width
    
       cdef int sizeout = np.ceil(rmax/dr)

       cdef np.ndarray[np.float64_t, ndim=1] radial = np.empty(sizeout,dtype=np.float64)
       cdef np.ndarray[np.float64_t, ndim=1] radialdatamean = np.empty(sizeout,dtype=np.float64)

       radial = np.arange(rmax/dr)*dr + np.float64(dr/2)
       nrad = len(radial)

       # Loop through the bins
       for irad from 0 <= irad < nrad:
         minrad = irad*dr
         maxrad = minrad + dr
         radialdatamean[irad] = data[(r>=minrad) * (r<maxrad)].mean()

       return radialdatamean

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef double _get_lengthscale(self, np.ndarray[np.float64_t, ndim=2] im, int lentype, double res, int nx, int ny):
      '''
      Calculates the integral lengthscale of 
      a tapered 2D grid of a Nx3 point cloud
      using 1 of 3 available methods

      Syntax
      ----------
      lengthscale = pysesa.lengthscale._get_lengthscale(points, lentype, res, nx, ny)

      Parameters
      ----------
      im : ndarray
   	   tapered 2D grid of Nx3 point cloud
      res : float
           spatial grid resolution to create a grid
      lentype : int
   	   lengthscale type:
           1, l<0.5
           2, l<1/e
           3, l<0
      nx : int
      	   size of im in x direction
      ny : int
      	   size of im in y direction

      Returns
      ----------
      self.lengthscale: float
   	   integral lengthscale

      '''
      cdef double ma
      cdef np.ndarray[np.complex128_t, ndim=2] ft1 = np.empty((nx,ny),dtype=np.complex128)
      cdef np.ndarray[np.float64_t, ndim=2] mag1 = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] auto = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] autoarray = np.empty((nx,ny),dtype=np.float64)

      # =========================================================
      # autocorrelation
      # 2D fourier transform of demeaned image
      if np.sum(np.isnan(im)) == np.prod(np.shape(im)):
         return np.nan

      if np.prod(np.shape(im))<2:
         return np.nan

      ft1 = np.fft.fft2(im)
      mag1 = pow(abs(ft1),2) # power spectrum
      # autocovariance as inverse fourier transform as zero-centred power spectrum
      auto = np.fft.fftshift(np.real(np.fft.ifft2(mag1)))
      autoarray = np.asarray(auto) # make 1d to find max
      ma = autoarray.max() # find max
      # get integral lengthscale
      if ma>0:
         auto = np.dot(auto,pow(ma,-1))
         h = self._radial_data(np.squeeze(auto))

         try:
            if lentype==1: #l<0.5
               return res*(2*np.pi)*np.abs(trapz(h[:np.where(h<0.5)[0][0]+1], np.arange(np.where(h<0.5)[0][0]+1) ))
            elif lentype==2: #l<1/e
               return res*(2*np.pi)*np.abs(trapz(h[:np.where(h<(1/np.exp(1)))[0][0]+1], np.arange(np.where(h<(1/np.exp(1)))[0][0]+1) ))
            else: #l<0
               return res*(2*np.pi)*np.abs(trapz(h[:np.where(h<0)[0][0]+1], np.arange(np.where(h<0)[0][0]+1) ))
         except:
            return np.nan
      else:
         return np.nan

#         try:
#            if lentype==1: #l<0.5
#               return res*((2*np.pi)*(np.where(h<0.5)[0][0]+1))
#            elif lentype==2: #l<1/e
#               return res*((2*np.pi)*(np.where(h<(1/np.exp(1)))[0][0]+1))
#            else: #l<0
#               return res*(np.where(h<0)[0][0]+1)
#         except:
#            return np.nan
#      else:
#         return np.nan

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] _Hanning2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      Returns a 2D Hanning (a weighted cosine) tapered 2D data

      Syntax
      ----------
      im = pysesa.lengthscale._Hanning2D(im)

      Parameters
      ----------
      im : ndarray
   	   2D grid of Nx3 point cloud

      Returns
      ----------
      im: ndarray
   	   tapered 2D grid 

      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.hanning(ny),np.hanning(nx)))
      return im*Wss

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] _Hamming2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      Returns a 2D Hamming (a weighted cosine) tapered 2D data

      Syntax
      ----------
      im = pysesa.lengthscale._Hamming2D(im)

      Parameters
      ----------
      im : ndarray
   	   2D grid of Nx3 point cloud

      Returns
      ----------
      im: ndarray
   	   tapered 2D grid 

      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.hamming(ny),np.hamming(nx)))
      return im*Wss

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] _Blackman2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      Returns a 2D Blackman (a summation of cosines) tapered 2D data

      Syntax
      ----------
      im = pysesa.lengthscale._Blackman2D(im)

      Parameters
      ----------
      im : ndarray
   	   2D grid of Nx3 point cloud

      Returns
      ----------
      im: ndarray
   	   tapered 2D grid 

      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.blackman(ny),np.blackman(nx)))
      return im*Wss

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] _Bartlett2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      Returns a 2D Bartlett (a triangular) tapered 2D data

      Syntax
      ----------
      im = pysesa.lengthscale._Bartlett2D(im)

      Parameters
      ----------
      im : ndarray
   	   2D grid of Nx3 point cloud

      Returns
      ----------
      im: ndarray
   	   tapered 2D grid 

      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.bartlett(ny),np.bartlett(nx)))
      return im*Wss

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] _taper_im(self, np.ndarray[np.float64_t, ndim=2] im, int taper):
      '''
      Returns the tapered 2D grid of 3D pointcloud for spectral analysis

      Syntax
      ----------
      im = pysesa.lengthscale._taper(im, taper)

      Parameters
      ----------
      im : ndarray
   	   2D grid of Nx3 point cloud
      taper : int
   	   flag for taper type:
           1, Hanning (Hann)
           2, Hamming
           3, Blackman
           4, Bartlett

      Returns
      ----------
      im: ndarray
   	   tapered 2D grid of 3D pointcloud

      '''
      im[np.isnan(im)] = np.nanmean(im.flatten())
      im[np.isinf(im)] = np.nanmean(im.flatten())
      im[im<=0] = np.nanmean(im.flatten())

      if taper==1:
         im = self._Hanning2D(im)
      elif taper==2:
         im = self._Hamming2D(im)
      elif taper==3:
         im = self._Blackman2D(im)
      else:
         im = self._Bartlett2D(im)
      return im

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef double getlengthscale(self):
      '''
      Returns the integral lengthscale

      Syntax
      ----------
      lengthscale = pysesa.lengthscale.getlengthscale()

      Parameters
      ----------
      self : instance
   	   pysesa.lengthscale instance

      Returns
      ----------
      self.lengthscale: float
   	   integral lengthscale

      '''
      return self.lengthscale


   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray[np.float64_t, ndim=2] getdata(self):
      '''
      Returns the tapered 2D grid of 3D pointcloud for spectral analysis

      Syntax
      ----------
      im = pysesa.lengthscale.getdata()

      Parameters
      ----------
      self : instance
   	   pysesa.lengthscale instance

      Returns 
      ----------
      self.data: ndarray
   	   tapered 2D grid of 3D pointcloud

      '''
      return self.data


