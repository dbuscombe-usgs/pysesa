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
##                _ __     
## _      _______(_) /____ 
##| | /| / / ___/ / __/ _ \
##| |/ |/ / /  / / /_/  __/
##|__/|__/_/  /_/\__/\___/ 
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
import numpy as np
cimport numpy as np
cimport cython

# =========================================================
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.nonecheck(False)
def txtwrite(str outfile, np.ndarray[np.float64_t, ndim=2] towrite, str header=None):
   '''
   Custom fast numpy array to comma-delimited ASCII txt file

   Syntax
   ----------
   () = pysesa.write.txtwrite(infile, towrite)

   Parameters
   ------------
   outfile : str
   	name of file to write to
   towrite : ndarray
   	ndarray containing Nx3 point cloud

   Other Parameters
   -----------------
   header : str, *optional* [default = None]
   	header string

   Returns
   ----------
   None

   '''
   with open(outfile, 'w') as f:

      np.savetxt(f, towrite[np.where(towrite[:,-1])[0],:], header = header, fmt=' '.join(['%8.6f,'] * np.shape(towrite)[1])[:-1]) 

