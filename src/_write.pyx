# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8

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
   with open(outfile, 'wb') as f:

      np.savetxt(f, towrite[np.where(towrite[:,-1])[0],:], header = header, fmt=' '.join(['%8.6f,'] * np.shape(towrite)[1])[:-1])
