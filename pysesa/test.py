## PySESA (Python program for Spatially Explicit Spectral Analysis)
## has been developed at the Grand Canyon Monitoring & Research Center,
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

import pysesa
import os
import shutil
import errno

__all__ = [
    'test',
    ]

def dircopy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else:
            print('Directory not copied. Error: %s' % e)

def test():
   """
   PySESA - a Python framework for Spatially Explicit Spectral Analysis

   PySESA is an open-source project dedicated to provide a generic Python framework
   for spatially explicit statistical analyses of point clouds and other geospatial data,
   in the spatial and frequency domains, for use in the geosciences

   The program is detailed in:
   Buscombe, D. (2016) "Computational considerations for spatially explicit spectral analysis of point clouds and geospatial data", 86, 92-108, 10.1016/j.cageo.2015.10.004.

   :Author:
       Daniel Buscombe
       Grand Canyon Monitoring and Research Center
       United States Geological Survey
       Flagstaff, AZ 86001
       dbuscombe@usgs.gov

   For more information visit http://dbuscombe-usgs.github.io/pysesa/

   :install:
       python setup.py install
       sudo python setup.py install

   :test:
       python -c "import pysesa; pysesa.test()"

   :license:
       GNU Lesser General Public License, Version 3
       (http://www.gnu.org/copyleft/lesser.html)

       This software is in the public domain because it contains materials that
       originally came from the United States Geological Survey, an agency of the
       United States Department of Interior. For more information,
       see the official USGS copyright policy at
       http://www.usgs.gov/visual-id/credit_usgs.html#copyright
       Any use of trade, product, or firm names is for descriptive purposes only
       and does not imply endorsement by the U.S. government.

   """
   # copy files over to somewhere read/writeable
   dircopy(pysesa.__path__[0], os.path.expanduser("~")+os.sep+'pysesa_test')
   shutil.copy(pysesa.__path__[0]+os.sep+'example_100000pts.xyz', os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_100000pts.xyz')

   # work on the 100,000 point data set
   infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_100000pts.xyz'

   out = 0.5 #m output grid
   detrend = 4 #ODR plane
   proctype = 4 #Processing spatial and spectral parameters (no smoothing)
   mxpts = 2048 # max pts per window
   res = 0.05 #cm internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # Hann taper
   prc_overlap = 100 # 300% overlap between successive windows
   minpts = 32 # min pts per window
   nchunks = 1 # split data into nchunks
   filt = 0 #no filter

   pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap, nchunks, filt)

   # work on the 2,000,000 point data set
   infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_2000000pts.xyz'
   nchunks = 2 # split data into nchunks
   filt = 1 #apply filter

   pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap, nchunks, filt)

if __name__ == '__main__':
   test()
