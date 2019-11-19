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
from scipy.integrate import trapz
from sklearn.neighbors import KDTree

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

   cdef object data, lengthscale, tree, rs1

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

      points = points - np.nanmin(points, axis=0)

      tree = KDTree(np.column_stack((points[:,0], points[:,1])))
      tmp = 1/tree.two_point_correlation(np.column_stack((points[:,0], points[:,1])), np.linspace(res,res*10,100))
      h = tmp/np.max(tmp)
      try:
        if lentype==1: #l<0.5
           out= res*(2*np.pi)*np.abs(trapz(h[:np.where(h<0.5)[0][0]+1], np.arange(np.where(h<0.5)[0][0]+1) ))
        elif lentype==2: #l<1/e
           out= res*(2*np.pi)*np.abs(trapz(h[:np.where(h<(1/np.exp(1)))[0][0]+1], np.arange(np.where(h<(1/np.exp(1)))[0][0]+1) ))
        else: #l<0
           if h.min()<0:
              out= res*(2*np.pi)*np.abs(trapz(h[:np.where(h<0)[0][0]+1], np.arange(np.where(h<0)[0][0]+1) ))
           else:
              out= res*(2*np.pi)*np.abs(trapz(h[:np.where(h<(1/np.exp(1)))[0][0]+1], np.arange(np.where(h<(1/np.exp(1)))[0][0]+1) ))
      except:
        out= np.nan

      self.lengthscale = out

      # pre-allocate some more arrays and get the ranges of x and y
      cdef float nnx = np.max(points[:,0]) - np.min(points[:,0])
      cdef float nny = np.max(points[:,1]) - np.min(points[:,1])
      cdef int lenx
      cdef int nx, ny, i

      # enfore square matrix in gridding
      if nny!=nnx:
         lenx = int(np.ceil(nnx/(res*(nnx/nny))))
      else:
         lenx = int(np.ceil(nnx/res))

      cdef np.ndarray[np.float64_t, ndim=2] grid_x = np.empty((lenx,lenx), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] grid_y = np.empty((lenx,lenx), dtype=np.float64)

      cdef np.ndarray[np.float64_t, ndim=2] im = np.empty((lenx,lenx), dtype=np.float64)

      _, inds = tree.query(np.column_stack((grid_x.flatten(), grid_y.flatten())), k = 1)
      im = points[:,2].flatten()[inds].reshape(np.shape(grid_x))
      im = self._taper_im(im, taper)
      self.data = im

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

      im[np.isnan(im)] = np.mean(im.flatten())
      im[np.isinf(im)] = np.mean(im.flatten())
      im[im<=0] = np.mean(im.flatten())

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
