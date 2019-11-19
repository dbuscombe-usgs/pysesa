## PySESA (Python program for Spatially Explicit Spectral Analysis)
## has been developed at the Grand Canyon Monitorinf & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
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
##                         __             __
##   _________  ___  _____/ /__________ _/ /
##  / ___/ __ \/ _ \/ ___/ __/ ___/ __ `/ /
## (__  ) /_/ /  __/ /__/ /_/ /  / /_/ / /
##/____/ .___/\___/\___/\__/_/   \__,_/_/
##    /_/

##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
##|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
##"""

# import libraries
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from scipy.integrate import trapz
from scipy.stats import linregress
import statsmodels.api as smapi

#from libc.math cimport sqrt,abs

# import PySESA libraries
import lengthscale as lengthscale_func

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')
import warnings
warnings.filterwarnings("ignore")

# =========================================================
cdef class spectral:
   '''
   Calculate spectral statistics of a Nx3 point cloud

   Syntax
   ----------
   data = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getdata()

   lengths = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengths()

   psdparams= pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getstats()

   lengthscale = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getlengthscale()

   moments = pysesa.spectral.spec(points, nbin, res, proctype, lentype, taper, method).getmoments()

   Parameters
   ------------
   points : ndarray
   	Nx3 point cloud

   Other Parameters
   -------------------
   nbin : int, *optional* [default = 20]
        number of bins for power spectral binning
   res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
   proctype : int, *optional* [default = 1, no spectral smoothing]
   	proctype type:
        1, no spectral smoothing

        2, spectrum smoothed with Gaussian

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
   ----------------------------------------
   self.data: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

        l = integral lengthscale

        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density

        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

        phi = effective slope (degrees)

   Returns [requested through .getpsdparams()]
   --------------------------------------------
   self.psdparams: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

   Returns [requested through .getlengths()]
   -------------------------------------------
   self.lengths: list
        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density


   Returns [requested through .getlengthscale()]
   ----------------------------------------------
   self.lengthscale: float
        l = integral lengthscale

   Returns [requested through .getmoments()]
   ------------------------------------------
   self.moments: list
        Z = zero-crossings per unit length

        E = extreme per unit length
        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

        phi = effective slope (degrees)

   '''
   cdef object data, x_space, k_space, r, moments, lengthscale, lengths, psdparams

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] points, int nbin=20, double res=0.05, int proctype=1, int lentype=1, int taper=1, str method='nearest'):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      data = pysesa.spectral(points, nbin, res, proctype, lentype, taper, method).getdata()

      lengths = pysesa.spectral(points, nbin, res, proctype, lentype, taper, method).getlengths()

      psdparams= pysesa.spectral(points, nbin, res, proctype, lentype, taper, method).getstats()

      lengthscale = pysesa.spectral(points, nbin, res, proctype, lentype, taper, method).getlengthscale()

      moments = pysesa.spectral(points, nbin, res, proctype, lentype, taper, method).getmoments()

      Parameters
      ------------
      points : ndarray
   	Nx3 point cloud

      Other Parameters
      ------------------
      nbin : int, *optional* [default = 20]
        number of bins for power spectral binning
      res : float, *optional* [default = 0.05]
        spatial grid resolution to create a grid
      proctype : int, *optional* [default = 1, no spectral smoothing]
   	proctype type:
        1, no spectral smoothing

        2, spectrum smoothed with Gaussian

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
      ---------------------------------------
      self.data: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

        l = integral lengthscale

        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density

        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

        phi = effective slope (degrees)

      Returns [requested through .getpsdparams()]
      ---------------------------------------------
      self.psdparams: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension


      Returns [requested through .getlengths()]
      -------------------------------------------
      self.lengths: list
        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density


      Returns [requested through .getlengthscale()]
      -----------------------------------------------
      self.lengthscale: float
        l = integral lengthscale

      Returns [requested through .getmoments()]
      -------------------------------------------
      self.moments: list
        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

      '''

      cdef float pi = 3.14159265

      cdef int nx, ny
      cdef double slope, intercept, r_value, p_value, std_err, l, rms1, rms2, wmax, Z, E, sigma, T0_1, T0_2, sw1, sw2, d, wmean
      cdef np.ndarray[np.float64_t, ndim=1] moment = np.empty(5,dtype=np.float64)

      cdef np.ndarray[np.float64_t, ndim=1] s
      cdef np.ndarray[np.float64_t, ndim=1] s_b
      cdef np.ndarray[np.float64_t, ndim=1] k_back
      cdef np.ndarray[np.float64_t, ndim=1] k

      #cdef double[::1] s
      cdef np.ndarray[np.float64_t, ndim=2] im

      points = points - np.nanmin(points, axis=0)

      r = lengthscale_func.lengthscale(points, res, lentype, taper, method)
      im = r.getdata()
      l = r.getlengthscale()
      self.lengthscale = l

      ny, nx= np.shape(im)

      # is all nans just return nans
      if np.sum(np.isnan(im)) == np.prod(np.shape(im)):

         self.data = [np.ones(24)*np.nan]
         return

      else:

         #if proctype==1: #no spectral smoothing

         #try:

         s, k  = self._psd(im, nbin, res, proctype) #k_back, s_b
         slope, intercept, r_value, p_value, std_err, d = self._psparams(k, s, res)
         self.psdparams = [slope, intercept, r_value, p_value, std_err, d]

         # interpolate background spectrum onto wavenumber array
         #s_b = np.interp(k,k_back,s_b)

         k = k[1:] #Nyquist frequency up
         s = s[1:]
         #s_b = s_b[1:]

         # wavelengths and rms amplitudes
         wmax, wmean, rms1, rms2 = self._wav_rms(k, s, res) #s_b
         self.lengths = [wmax, wmean, rms1, rms2]

         # moments and moment parameters
         Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = self._moments(k, s, res)
         self.moments = [Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4]

         # concetanate and add effective slope
         #self.data = self.psdparams + [self.lengthscale] + self.lengths + self.moments + [np.arctan(sigma/l)/(np.pi/180)]
         self.data = self.psdparams + [self.lengthscale] + self.lengths + self.moments + [np.arctan(sigma/l)/(pi/180)]

         #except:

        #    self.data = [np.ones(24)*np.nan]

         return

#         elif proctype==2: # with spectral smoothing

#            try:

#               s, s_b, k, k_back = self._psd(im, nx, ny, nbin, res, proctype)
#               slope, intercept, r_value, p_value, std_err, d = self._psparams(k_back, s_b, res)
#               self.psdparams = [slope, intercept, r_value, p_value, std_err, d]
#
#               # interpolate background spectrum onto wavenumber array
#               s_b = np.interp(k,k_back,s_b)

#               k = k[1:]
#               s = s[1:]
#               s_b = s_b[1:]

#               # wavelengths and rms amplitudes
#               wmax, wmean, rms1, rms2 = self._wav_rms(k, s, s_b, res)
#               self.lengths = [wmax, wmean, rms1, rms2]
#
#               # moments and moment parameters
#               Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = self._moments(k, s_b, res)
#               self.moments = [Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4]

#               # concetanate and add effective slope
#               #self.data = self.psdparams + [self.lengthscale] + self.lengths + self.moments + [np.arctan(sigma/l)/(np.pi/180)]
#               self.data = self.psdparams + [self.lengthscale] + self.lengths + self.moments + [np.arctan(sigma/l)/(pi/180)]

#            except:

#               self.data = [np.ones(24)*np.nan]

#            return

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray _movingaverage(self, np.ndarray[np.float64_t, ndim=1] interval, int window_size):
      return np.convolve(interval, np.ones(int(window_size))/float(window_size), 'same')

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _moments(self, np.ndarray[np.float64_t, ndim=1] k, np.ndarray[np.float64_t, ndim=1] s_b, double res):
      '''
      Return moments and moment parameters

      Syntax
      ----------
      Z, E, sigma, T0_1, T0_2, sw1, sw2, m0, m1, m2, m3, m4 = pysesa_spectral.spec._moments(k, s_b, res)

      Parameters
      ------------
      k : ndarray
   	   N x 1 wavenumbers
      s_b : ndarray
   	   N x 1 array of background power spectral density
      res : double
   	   spatial resolution of grid

      Returns
      ----------
      data: tuple of floats
        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

      '''
      cdef double Z, E, sigma, T0_1, T0_2, sw1, sw2
      cdef np.ndarray[np.float64_t, ndim=1] moment = np.empty(5,dtype=np.float64)
      cdef int i

      # get moments of spectrum
      for i from 0 <= i < 5:
      #for i in xrange(0,5):
         moment[i] = np.abs(trapz((k)**i,s_b)) #,np.median(np.gradient(k))))
         #moment[i] = abs(trapz((k)**i,s_b)) #,np.median(np.gradient(k))))

      # s is actually 10**length^4
      moment = np.sqrt((10**4)*moment)
      #moment = sqrt((10**4)*moment)

      # zero crossings per second
      Z = 2*(moment[2]/moment[0])**0.5
      #Z = 2*sqrt(moment[2]/moment[0])

      # extrema per second
      E = 2*(moment[4]/moment[2])**0.5
      #E = 2*sqrt(moment[4]/moment[2])

      #rms
      sigma = (moment[2]/moment[0])

      # average period m0/m1
      T0_1 = (moment[0]/moment[1])
      #average period (m0/m2)^0.5
      T0_2 = T0_1**0.5
      # spectral width parameter
      sw1 = (moment[0]*moment[2]/moment[1]**2-1)**0.5
      # spectral width paramenter
      sw2 = np.abs(1 - moment[2]**2/(moment[0]*moment[4]))**0.5
      #sw2 = abs(1 - moment[2]**2/(moment[0]*moment[4]))**0.5

      return [Z, E, sigma, T0_1, T0_2, sw1, sw2] + moment.tolist()

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _wav_rms(self, np.ndarray[np.float64_t, ndim=1] k, np.ndarray[np.float64_t, ndim=1] s, double res): #np.ndarray[np.float64_t, ndim=1] s_b,
      '''
      Return max and mean wavelengths and rms amplitudes

      Syntax
      ----------
      wmax, wmean, rms1, rms2 = pysesa_spectral.spec._wav_rms(k, s, res) #s_b,

      Parameters
      ----------
      k : ndarray
   	   N x 1 wavenumbers
      s : ndarray
   	   N x 1 array of power spectral density
      s_b : ndarray
   	   N x 1 array of background power spectral density
      res : double
   	   spatial resolution of grid

      Returns
      ----------
      data: tuple of floats
        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density

      '''
      cdef double wmax, wmean, rms1, rms2

      #cdef float pi = 3.14159265
      cdef np.ndarray[np.float64_t, ndim=1] ssb = np.empty(len(s),dtype=np.float64)

      s_b = self._movingaverage(s, len(s)/4)

      ssb = self._movingaverage(s/s_b,len(s)/10) #s_b

      # get peak wavelength
      wmax = res*(2*np.pi)/k[np.argmax(ssb)]
      #wmax = res*(2*np.pi)/k[np.argmax(np.abs(s/s_b))]
      #wmax = res*(2*pi)/k[np.argmax(abs(s/s_b))]

      # mean wavelength
      wmean = res*(2*np.pi)/np.abs(trapz(s/s_b, k))
      #wmean = res*(2*pi)/abs(trapz(s/s_b, k))

      # get rms amplitudes
      rms1 = np.sqrt(np.abs(trapz(s, k)))/res
      rms2 = np.sqrt(np.abs(trapz(s/ssb, k)))/res
      #rms2 = np.sqrt(np.abs(trapz(s_b, k)))/res
      #rms1 = sqrt(abs(trapz(s, k)))/res
      #rms2 = sqrt(abs(trapz(s_b, k)))/res

      return [wmax, wmean, rms1, rms2]

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _psd(self, np.ndarray[np.float64_t, ndim=2] im, int nbin, double res, int proctype):
      '''
      Return spectrum, wavenumber index, and background spectrum

      Syntax
      ----------
      s, s_b, k, k_back = pysesa_spectral.spec._psd(im, nx, ny, nbin, res, proctype)

      Parameters
      ----------
      im : ndarray
   	   2D array
      nx : int
   	   size of im in x dimension
      ny : int
   	   size of im in y dimension
      nbin : int
   	   number of bins for spectral binning
      res : double
   	   spatial resolution of grid
      proctype : int
   	   1 = no smoothing, 2 = smoothing

      Returns
      ----------
      data: tuple of ndarrays
   	s = power spectral density

        s_b = background power spectral density

        k = wavenumber index of s

        k_back = wavenumber index of s_b

      '''

      #from nifty import rg_space, field, about
      #about.warnings='OFF'
      im = im - np.nanmean(im.flatten()) #detrend(im)

      cdef int nxh
      ny, nx= np.shape(im)
      nxh = int(nx/2)+1

      #cdef int dfx , dfy
      cdef np.float eps
      cdef int xc, yc
      cdef double dfx, dfy, ma
      cdef np.ndarray[np.float64_t, ndim=1] Wss = np.empty(nx,dtype=np.float64)
      cdef np.ndarray[np.complex128_t, ndim=2] M = np.empty((nx,ny),dtype=np.complex128)
      cdef np.ndarray[np.int32_t, ndim=2] rows = np.empty((nx,ny),dtype=np.int32)
      cdef np.ndarray[np.int32_t, ndim=2] cols = np.empty((nx,ny),dtype=np.int32)
      cdef np.ndarray[np.float64_t, ndim=2] fm = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.complex128_t, ndim=2] Pm = np.empty((nx,ny),dtype=np.complex128)
      cdef np.ndarray[np.complex128_t, ndim=2] Pmc = np.empty((nx,nxh),dtype=np.complex128)
      cdef np.ndarray[np.float64_t, ndim=2] fv = np.empty((nx,nxh),dtype=np.float64)

      cdef np.ndarray[np.complex128_t, ndim=2] a #= np.empty((nx*nxh,2),dtype=np.float128)
      cdef np.ndarray[np.float64_t, ndim=1] Pv #= np.empty(nx*nxh,dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=1] fv2 #= np.empty(nx*nxh,dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] B

      cdef np.ndarray[np.float64_t, ndim=1] s
      #cdef np.ndarray[np.float64_t, ndim=1] s_b
      #cdef np.ndarray[np.float64_t, ndim=1] k_back
      cdef np.ndarray[np.float64_t, ndim=1] k
      #cdef object x_space, k_space

      #s, s_b, k, k_back = self._psd(im, nx, ny, nbin, res, proctype)

      im[np.isnan(im)] = np.nanmean(im.flatten())

      eps = 2.22e-16
      dfx =1/np.asarray(nx).astype(float)
      dfy = 1/np.asarray(ny).astype(float)

      im, Wss = self._Hann2D(im) # window the matrix with an elliptical Hann (raised cosine) window

      # 2D fourier transform of demeaned image
      M= np.fft.fftshift(np.fft.fft2(im - im.mean())) # sp.detrend(im)))# im-imMean))
      # zero out the DC component
      M[int(ny/2) + 1, int(nx/2) + 1] = 0

      xc = nx/2+1; yc = ny/2+1
      cols, rows = np.meshgrid(np.r_[:nx],np.r_[:ny])
      # frequency matrix
      fm = np.sqrt((dfy*(rows-yc))**2 + (dfx*(cols-xc))**2)
      # Calculate the DFT periodogram, which has units of amplitude^2
      # Dividing by Wss*twopower^2 corrects for the reduction of
      # amplitude by the windowing function
      Pm = M * np.conj(M) / (nx * ny * Wss)

      # get 1D spectrum
      # Create sorted, non-redundant vectors of frequency and power
      Pmc = Pm[:,:nxh]
      fv = fm[:,:nxh]
      fv[ np.asarray(np.r_[(yc+1):ny],'int') , : ] = -1 # % This half-column is redundant.
      a = np.vstack( (fv.flatten(), Pmc.flatten()) ).T
      a = a[a[:,0].argsort(),]	#Sort rows (by first row)
      a = a[a[:,0]>0,:]
      Pv = np.real(2*a[:,1])
      fv2 = np.real(a[:,0])

      # remove tiny powers less than machine precision
      index = np.where(Pv>eps)
      s = Pv[index]
      k = fv2[index]*1/res

    # # interpolate background spectrum onto wavenumber array
    # s_b = np.interp(k,k_back,s_b)
    #
    # k = k[2:] #Nyquist frequency up
    # s = s[2:]
    # s_b = s_b[2:]

      # # set up field
      # if nx%2==0:
      #    x_space = rg_space((nx,ny), naxes=2)
      # else:
      #    x_space = rg_space((nx-1,ny-1), naxes=2)
      # k_space = x_space.get_codomain()
      #
      # if proctype==2:
      #    # and get the smoothed power spectrum
      #    if nx%2==0:
      #       s = field(x_space, target=k_space, val=im).power(smooth=1)
      #    else:
      #       s = field(x_space, target=k_space, val=im[:-1,:-1].T).power(smooth=1)
      # else:
      #    # and get the power spectrum
      #    if nx%2==0:
      #       s = field(x_space, target=k_space, val=im).power(smooth=0)
      #    else:
      #       s = field(x_space, target=k_space, val=im[:-1,:-1].T).power(smooth=0)
      #
      # k = k_space.power_indices["kindex"]*res
      # # generate synthetic field and get binned power spectrum
      # s_b = field(x_space, target=k_space, random="syn", spec=s).power(nbin=nbin)
      # k_back = k_space.power_indices["kindex"]*res

      return [s, k] #s_b, k_back

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   def _Hann2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      return a 2D Hann taper
      '''
      cdef int nx, ny, a, b

      ny, nx= np.shape(im)

      #cdef np.ndarray[np.int64_t, ndim=2] X = np.empty((nx,ny),dtype=np.int64)
      #cdef np.ndarray[np.int64_t, ndim=2] Y = np.empty((nx,ny),dtype=np.int64)
      cdef np.ndarray[np.float64_t, ndim=2] theta = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] r = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] rprime = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] hancoeff = np.empty((nx,ny),dtype=np.float64)

      cdef np.ndarray[np.float64_t, ndim=2] res = np.empty((nx,ny),dtype=np.float64)

      a = nx/2+1
      b = ny/2+1
      X, Y = np.meshgrid(np.r_[:nx],np.r_[:ny])
      theta = np.asarray(X==a,'int') * (np.pi/2) + np.asarray(X!=a,'int') * np.arctan2( (Y-b),(X-a) )
      # angular polar coordinate
      r = np.sqrt((Y-b)**2 + (X-a)**2) # radial polar coordinate
      rprime = np.sqrt((a**2)*(b**2)*(b**2*(np.cos(theta))**2 + a**2*(np.sin(theta))**2)**(-1))
      # 'radius' of ellipse for this theta
      hanncoeff = (r < rprime)*(0.5*(1 + np.cos(np.pi*r/rprime)))

      res= np.asarray(im*hanncoeff,np.float64)
      return res, sum(hanncoeff**2)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list _psparams(self, np.ndarray[np.float64_t, ndim=1] k, np.ndarray[np.float64_t, ndim=1] s, double res):
      '''
      Return slope, intercept, r_value, p_value, std_err, fractal dimension

      Syntax
      ----------
      slope, intercept, r_value, p_value, std_err, d = pysesa_spectral.spec._psdparams(k, s, res)

      Parameters
      ------------
      k_back : ndarray
   	   N x 1 wavenumbers
      s_b : ndarray
   	   N x 1 array of power spectral density
      res : double
   	   spatial resolution of grid

      Returns
      ----------
      data: tuple
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

      '''
      cdef double slope, intercept, r_value, p_value, std_err, d

      slope, intercept, r_value, p_value, std_err = self._do_linreg(np.c_[np.log10(k),np.log10(s)])
      std_err = std_err*res

      ## give intercept square units
      intercept = (10**intercept)**(1/(2-slope))

      if intercept>1:
         intercept=1
         slope=0

      if slope>0:
         slope=0

      # fractal dimension
      d = (8+slope)/2

      return [slope, intercept, r_value, p_value, std_err, d]

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef tuple _do_linreg(self, np.ndarray[np.float64_t, ndim=2] B):
      '''
      Do a robust linear regression

      Syntax
      ----------
      slope, intercept, r_value, p_value, std_err = pysesa_spectral.spec._do_linreg(B)

      Parameters
      ------------
      B : ndarray
   	   N x 2 array of wavenumber and power spectral density

      Returns
      ----------
      data: tuple
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

      '''

      cdef float slope, intercept, r_value, p_value, std_err
      cdef object regression

      # remove any rows with nans and infs
      B = B[np.where(np.logical_not(np.any(np.isinf(B),axis=1)))[0],:]
      B = B[np.where(np.logical_not(np.any(np.isnan(B),axis=1)))[0],:]

      # get OLS regression
      slope, intercept, r_value, p_value, std_err = linregress(B[:,0], B[:,1])
      # try to RLM regression for slope and intercept
      try:
         regression = smapi.RLM(B[:,1], smapi.add_constant(B[:,0], prepend=True) ).fit()
         intercept = regression.params[0]
         slope = regression.params[1]
      except:
         pass
      return (slope, intercept, r_value, p_value, std_err)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      data = pysesa.spectral.getdata()

      Parameters
      ------------
      self : instance
   	   pysesa.spectral instance

      Returns [requested through .getdata()]
      -----------------------------------------
      self.data: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        l = integral lengthscale

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density

        wmax = peak wavelength

        wmean = mean wavelength

        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        d = fractal dimension

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

        phi = effective slope (degrees)
      '''
      return self.data

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getmoments(self):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      moments = pysesa.spectral.getmoments()

      Parameters
      ------------
      self : instance
   	   pysesa.spectral instance

      Returns [requested through .getmoments()]
      ------------------------------------------
      self.moments: list
        Z = zero-crossings per unit length

        E = extreme per unit length

        sigma = RMS amplitude

        T0_1 = average spatial period (m_0/m_1)

        T0_2 = average spatial period (m_0/m_2)^0.5

        sw1 = spectral width

        sw2 = spectral width (normalised radius of gyration)

        m0 = zeroth moment of spectrum

        m1 = first moment of spectrum

        m2 = second moment of spectrum

        m3 = third moment of spectrum

        m4 = fourth moment of spectrum

      '''
      return self.moments

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef double getlengthscale(self):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      lengthscale = pysesa.spectral.getlengthscale()

      Parameters
      ------------
      self : instance
   	   pysesa.spectral instance

      Returns [requested through .getlengthscale()]
      -----------------------------------------------
      self.lengthscale: float
        l = integral lengthscale
      '''
      return self.lengthscale

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getlengths(self):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      lengths = pysesa.spectral.getlengths()

      Parameters
      ----------
      self : instance
   	   pysesa.spectral instance

      Returns [requested through .getlengths()]
      ------------------------------------------
      self.lengths: list
        wmax = peak wavelength

        wmean = mean wavelength

        rms1 = RMS amplitude from power spectral density

        rms2 = RMS amplitude from bin averaged power spectral density

      '''
      return self.lengths

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getpsdparams(self):
      '''
      Calculate spectral statistics of a Nx3 point cloud

      Syntax
      ----------
      psdparams= pysesa.spectral.getstats()

      Parameters
      ------------
      self : instance
   	   pysesa.spectral.spec instance

      Returns [requested through .getpsdparams()]
      ---------------------------------------------
      self.psdparams: list
   	slope = slope of regression line through log-log 1D power spectral density

        intercept = intercept of regression line through log-log 1D power spectral density

        r_value = correlation of regression through log-log 1D power spectral density

        p_value = probability that slope of regression through log-log 1D power spectral density is not zero

        std_err = standard error of regression through log-log 1D power spectral density

        d = fractal dimension

      '''
      return self.psdparams
