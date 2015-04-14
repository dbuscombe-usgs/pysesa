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
##                        __
##   ________  ____ _____/ /
##  / ___/ _ \/ __ `/ __  / 
## / /  /  __/ /_/ / /_/ /  
##/_/   \___/\__,_/\__,_/   


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
import numpy as np
cimport numpy as np
cimport cython

# =========================================================
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray txtread(str infile):
   '''
   Custom fast (up to 3.5x faster than numpy's genfromtxt) txt file to numpy array
   accepts comma, tab or space delimited files of 3 columns: x, y, and amplitude


   Syntax
   ----------
   pts = pysesa_read.txtread(infile)

   Parameters
   ------------
   infile : str
   	3-column ASCII file containing Nx3 point cloud

   Returns
   ----------
   data: ndarray
   	Nx3 point cloud, 32 bit precision

   '''
   f = open(infile, 'rb'); data = f.read(); f.close()
   data = data.splitlines()
   
   cdef np.ndarray[np.float32_t, ndim=2] out = np.empty((len(data),3),dtype=np.float32)
   
   if len(data[0].split(','))>1: # the file is comma delimited
      out = np.array([x.split(',',2)[0:3] for x in data], dtype='float32')
      return out
   else:
      try: # space delimited
         out = np.array([x.split(' ',2)[0:3] for x in data], dtype='float32')
         return out
      except: #tab delimited
         out = np.array([x.split('\t',2)[0:3] for x in data], dtype='float32')
         return out
         
