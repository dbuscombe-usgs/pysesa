"""                       
   _________  ___  _____
  / ___/ __ \/ _ \/ ___/
 (__  ) /_/ /  __/ /__  
/____/ .___/\___/\___/  
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

from scipy.integrate import trapz
from scipy.stats import linregress
import statsmodels.api as smapi

import warnings
warnings.filterwarnings("ignore")

# =========================================================
cdef class spec:
   """
   Returns an instance. All spectral parameters requested through .getdata()
   """
   cdef object data, x_space, k_space

   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] im, int nbin, double res, int proctype=1, int lentype=1, int taper=1): 
      '''
      get spectral parameters from point cloud
      '''

      cdef int nx, ny, i
      im = im - np.nanmean(im)

      ny, nx= np.shape(im)

      cdef double ma
      cdef np.ndarray[np.complex128_t, ndim=2] ft1 = np.empty((nx,ny),dtype=np.complex128)
      cdef np.ndarray[np.float64_t, ndim=2] mag1 = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] auto = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] autoarray = np.empty((nx,ny),dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((nx,ny),dtype=np.float64)

      cdef double slope, intercept, r_value, p_value, std_err, out, rms1, rms2, w, Z, E, sigma, T0_1, T0_2, sw1, sw2
      cdef np.ndarray[np.float64_t, ndim=1] moment = np.empty(5,dtype=np.float64)

      cdef np.ndarray[np.float64_t, ndim=1] s 
      cdef np.ndarray[np.float64_t, ndim=1] s_b 
      cdef np.ndarray[np.float64_t, ndim=1] k_back 

      # is all nans just return nans
      if np.sum(np.isnan(im)) == np.prod(np.shape(im)):

         slope = np.nan
         intercept = np.nan
         r_value = np.nan
         p_value = np.nan
         std_err = np.nan
         out = np.nan
         rms1 = np.nan
         rms2 = np.nan
         Z = np.nan
         E = np.nan										
         sigma = np.nan
         T0_1 = np.nan
         T0_2 = np.nan
         sw1 = np.nan
         sw2 = np.nan
         w = np.nan
         moment = moment * np.nan

         self.data = [slope, intercept, r_value, p_value, std_err, out, rms1, rms2, w, Z, E, sigma, T0_1, T0_2, sw1, sw2] + moment.tolist()
         return

      else:

         im[np.isnan(im)] = np.nanmean(im.flatten())
         im[np.isinf(im)] = np.nanmean(im.flatten())

         if taper==1:
            im, Wss = self._Hanning2D(im)
         elif taper==2:
            im, Wss = self._Hamming2D(im)
         elif taper==3:
            im, Wss = self._Blackman(im)
         else:
            im, Wss = self._Bartlett(im)

         # =========================================================
         # autocorrelation
         # 2D fourier transform of demeaned image
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
               if lentype==1:
                  out = res*((2*np.pi)*(np.where(h<0.5)[0][0]+1))
               else:
                  out = res*(np.where(h<0)[0][0]+1)
            except:
               out = np.nan
         else:
            out = np.nan

         if proctype==1: #no spectral smoothing

            from nifty import rg_space, field
            try:

               # set up field and get the power spectrum
               if nx%2==0:
                  x_space = rg_space((nx,ny), naxes=2)
               else:
                  x_space = rg_space((nx-1,ny-1), naxes=2)
               k_space = x_space.get_codomain()

               if nx%2==0:
                  s = field(x_space, target=k_space, val=im).power(smooth=0) # actual power spectrum
               else:
                  s = field(x_space, target=k_space, val=im[:-1,:-1]).power(smooth=0) # actual power spectrum

               k = k_space.power_indices["kindex"]*res #(res**2)
               # generate synthetic field and get binned power spectrum
               s_b = field(x_space, target=k_space, random="syn", spec=s).power(nbin=nbin) #field.power(a) 
               k_back = k_space.power_indices["kindex"]*res #(res**2)

               slope, intercept, r_value, p_value, std_err = self._do_linreg(np.c_[np.log10(k_back),np.log10(s_b)])
               std_err = std_err*res

               # interpolate background spectrum onto wavenumber array
               s_b = np.interp(k,k_back,s_b)

               # get peak wavelength
               w = res*((2*np.pi)/k[np.argmax(np.abs(s-s_b))]) #res*((2*np.pi)/k[np.argmax(s/s_b)])

               # get rms amplitudes
               rms1 = np.sqrt(np.abs(trapz(s, k)))/res
               rms2 = np.sqrt(np.abs(trapz(s_b, k)))/res

               # get moments of spectrum
               for i from 0 <= i < 5:
               #for i in xrange(0,5):
                  moment[i] = np.abs(trapz((k)**i,s_b,np.median(np.gradient(k))))

               Z = res*(2*np.sqrt(moment[2]/moment[0])) # zero crossings per second
               E = 2*np.sqrt(moment[4]/moment[2]) # extrema per second										
               sigma = (moment[2]/moment[0])*res #stdev

               T0_1 = (moment[0]/moment[1]) #/res    # average period m0/m1
               T0_2 = T0_1**0.5  #average period (m0/m2)^0.5
               sw1 = np.abs(moment[0]*moment[2]/moment[1]**2-1)**0.5 # spectral width parameter
               sw2 = np.abs(1 - moment[2]**2/(moment[0]*moment[4]))**0.5 # spectral width paramenter

            except:
                     
               slope = np.nan
               intercept = np.nan
               r_value = np.nan
               p_value = np.nan
               std_err = np.nan
               out = np.nan
               rms1 = np.nan
               rms2 = np.nan
               w = np.nan
               Z = np.nan
               E = np.nan										
               sigma = np.nan
               T0_1 = np.nan
               T0_2 = np.nan
               sw1 = np.nan
               sw2 = np.nan

            self.data = [slope, intercept, r_value, p_value, std_err, out, rms1, rms2, w, Z, E, sigma, T0_1, T0_2, sw1, sw2] + moment.tolist()
            return

         elif proctype==2: # with spectral smoothing

            from nifty import rg_space, field
            try:

               # set up field and get the smoothed power spectrum
               if nx%2==0:
                  x_space = rg_space((nx,ny), naxes=2)
               else:
                  x_space = rg_space((nx-1,ny-1), naxes=2)
               k_space = x_space.get_codomain()

               if nx%2==0:
                  s = field(x_space, target=k_space, val=im).power(smooth=1) # actual power spectrum
               else:
                  s = field(x_space, target=k_space, val=im[:-1,:-1]).power(smooth=1) # actual power spectrum

               k = k_space.power_indices["kindex"]*res #(res**2)
               # generate synthetic field and get binned power spectrum
               s_b = field(x_space, target=k_space, random="syn", spec=s).power(nbin=nbin) #field.power(a) 
               k_back = k_space.power_indices["kindex"]*res #(res**2)

               slope, intercept, r_value, p_value, std_err = self._do_linreg(np.c_[np.log10(k_back),np.log10(s_b)])
               std_err = std_err*res

               # interpolate background spectrum onto wavenumber array
               s_b = np.interp(k,k_back,s_b)

               # get peak wavelength
               w = res*((2*np.pi)/k[np.argmax(np.abs(s-s_b))]) #res*((2*np.pi)/k[np.argmax(s/s_b)])

               # get rms amplitudes
               rms1 = np.sqrt(np.abs(trapz(s, k)))/res
               rms2 = np.sqrt(np.abs(trapz(s_b, k)))/res

               # get moments of spectrum
               for i from 0 <= i < 5:
               #for i in xrange(0,5):
                  moment[i] = np.abs(trapz((k)**i,s_b,np.median(np.gradient(k))))# np.gradient(k)))

               Z = res*(2*np.sqrt(moment[2]/moment[0])) # zero crossings per second
               E = 2*np.sqrt(moment[4]/moment[2]) # extrema per second										
               sigma = (moment[2]/moment[0])*res #stdev

               T0_1 = (moment[0]/moment[1]) #/res    # average period m0/m1
               T0_2 = T0_1**0.5  #average period (m0/m2)^0.5
               sw1 = np.abs(moment[0]*moment[2]/moment[1]**2-1)**0.5 # spectral width parameter
               sw2 = np.abs(1 - moment[2]**2/(moment[0]*moment[4]))**0.5 # spectral width paramenter

            except:
                     
               slope = np.nan
               intercept = np.nan
               r_value = np.nan
               p_value = np.nan
               std_err = np.nan
               out = np.nan
               rms1 = np.nan
               rms2 = np.nan
               w = np.nan
               Z = np.nan
               E = np.nan										
               sigma = np.nan
               T0_1 = np.nan
               T0_2 = np.nan
               sw1 = np.nan
               sw2 = np.nan

            self.data = [slope, intercept, r_value, p_value, std_err, out, rms1, rms2, w, Z, E, sigma, T0_1, T0_2, sw1, sw2] + moment.tolist()
            return

         else: #proctype==3

            self.data = [out]
            return

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef np.ndarray _radial_data(self, np.ndarray[np.float64_t, ndim=2] data, int annulus_width=1):
       '''
       efficient radial average of matrix
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
   cpdef tuple _Hanning2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      return a 2D Hanning (a weighted cosine) taper
      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.hanning(ny),np.hanning(nx)))
      return (im*Wss, Wss)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef tuple _Hamming2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      return a 2D Hamming (a weighted cosine) taper
      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.hamming(ny),np.hamming(nx)))
      return (im*Wss, Wss)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef tuple _Blackman2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      return a 2D Blackman (a summation of cosines) taper
      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.blackman(ny),np.blackman(nx)))
      return (im*Wss, Wss)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef tuple _Bartlett2D(self, np.ndarray[np.float64_t, ndim=2] im):
      '''
      return a 2D Bartlett (a triangular) taper
      '''
      cdef int nx, ny
      ny, nx= np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] Wss = np.empty((ny, nx),dtype=np.float64)

      Wss = np.sqrt(np.outer(np.bartlett(ny),np.bartlett(nx)))
      return (im*Wss, Wss)

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef list getdata(self):
      return self.data

   # =========================================================
   @cython.boundscheck(False)
   @cython.cdivision(True)
   @cython.wraparound(False)
   @cython.nonecheck(False)
   cpdef tuple _do_linreg(self, np.ndarray[np.float64_t, ndim=2] B):
      '''
      do a robust linear regression
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

