# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8

import pysesa_main
import os

def test():
   """
   PySESA - a Python framework for Spatially Explicit Spectral Analysis

   PySESA is an open-source project dedicated to provide a generic Python framework
   for spatially explicit statistical analyses of point clouds and other geospatial data,
   in the spatial and frequency domains, for use in the geosciences

   The program is detailed in:
   Buscombe, D. (2016) "spatially explicit spectral analysis of point clouds and geospatial data", Computers and Geosciences 86, 92-108, 10.1016/j.cageo.2015.10.004.

   :Author:
       Daniel Buscombe
       United States Geological Survey
       Flagstaff, AZ 86001
       daniel.buscombe@nau.edu

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

   # work on the 100,000 point data set
   infile = 'data'+os.sep+'example_100000pts.xyz'

   out = 1 #m output grid
   detrend_mode = 4 #ODR plane\
   #detrend_mode = 1 #global mean
   mxpts = 2048 # max pts per window
   res = 0.01 #cm internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # Hann taper
   #prc_overlap = 100 # 100% overlap between successive windows
   prc_overlap = 100
   minpts = 16 # min pts per window
   nchunks = 1 # split data into nchunks (for large datasets)
   filt = 1 #apply filter

   proctype = 2 # spatial only

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                           lentype, minpts, taper, prc_overlap, nchunks, filt)

   # # work on the 2,000,000 point data set
   infile = 'data'+os.sep+'example_2000000pts.xyz'
   nchunks = 2 # split data into nchunks
   filt = 1 #apply filter

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                      lentype, minpts, taper, prc_overlap, nchunks, filt)

if __name__ == '__main__':
   test()
