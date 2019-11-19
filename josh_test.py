# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8

import pysesa_main
import os

def test():

   #################### file 1
   infile = 'josh_test'+os.sep+'Indoor2_Plot.las'

   out = 0.02 #m output grid
   detrend_mode = 4 #ODR plane\
   #detrend_mode = 1 #global mean
   mxpts = 1024 #2048 # max pts per window
   res = 0.01 #m internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # Hann taper
   prc_overlap = 100
   minpts = 8 #16 # min pts per window
   nchunks = 1 # split data into nchunks (for large datasets)
   filt = 1 #apply filter

   #################### file 2
   proctype = 2 # spatial only

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                           lentype, minpts, taper, prc_overlap, nchunks, filt)

   proctype = 1 # spectral only

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                           lentype, minpts, taper, prc_overlap, nchunks, filt)


   #################### file 2
   infile = 'josh_test'+os.sep+'P204sandy_2015_Focus_Nr_Undis.las'

   proctype = 2 # spatial only

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                           lentype, minpts, taper, prc_overlap, nchunks, filt)

   proctype = 1 # spectral only

   pysesa_main.process_all(infile, out, detrend_mode, proctype, mxpts, res, nbin,
                           lentype, minpts, taper, prc_overlap, nchunks, filt)


if __name__ == '__main__':
   test()
