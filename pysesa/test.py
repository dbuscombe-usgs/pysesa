"""
    ____  __  __________  _________ _
   / __ \/ / / / ___/ _ \/ ___/ __ `/
  / /_/ / /_/ (__  )  __(__  ) /_/ / 
 / .___/\__, /____/\___/____/\__,_/  
/_/    /____/   

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

#python -c "import pysesa; pysesa.test.dotest()"

import pysesa
import os
import shutil
import errno
 
def dircopy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else:
            print('Directory not copied. Error: %s' % e)

__all__ = [
    'dotest',
    ]

def dotest():

   # copy files over to somewhere read/writeable
   dircopy(pysesa.__path__[0], os.path.expanduser("~")+os.sep+'pysesa_test')
   shutil.copy(pysesa.__path__[0]+os.sep+'x_y_z_25cm.xyz', os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'test.DAT')

   # general settings   
   infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'x_y_z_25cm.xyz' 

   out = 1 #m output grid
   order = 3 #ODR plane
   proctype = 1 #Processing focal stats, lengthscale and spectral parameters (no smoothing)
   mxpts = 256 # max pts per window
   res = 0.05 #cm internal grid resolution
   nbin = 20 #number of bins for spectral binning
   lentype = 1 # l<0.5
   taper = 1 # yes do taper
   prc_overlap = 0 # no overlap between successive windows

   pysesa.pysesa(infile, out, order, proctype, mxpts, res, nbin, lentype, taper, prc_overlap)
   
if __name__ == '__main__':
   dotest()


