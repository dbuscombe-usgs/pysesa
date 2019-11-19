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

##https://github.com/grantbrown/laspy
from laspy.file import File

# ==== functions for rescaling las/laz data
def ascol( arr ):
    '''
    reshapes row matrix to be a column matrix (N,1).
    '''
    if len( arr.shape ) == 1: arr = arr.reshape( ( arr.shape[0], 1 ) )
    return arr

def scaled_x_dimension(las_file):
    x_dimension = las_file.X
    scale = las_file.header.scale[0]
    offset = las_file.header.offset[0]
    return(x_dimension*scale + offset)

def scaled_y_dimension(las_file):
    y_dimension = las_file.Y
    scale = las_file.header.scale[1]
    offset = las_file.header.offset[1]
    return(y_dimension*scale + offset)

def scaled_z_dimension(las_file):
    z_dimension = las_file.Z
    scale = las_file.header.scale[2]
    offset = las_file.header.offset[2]
    return(z_dimension*scale + offset)

# =========================================================
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray lasread(str infile):
   '''
   Laz/las file x,y,z to numpy array (just for elevation point clouds)


   Syntax
   ----------
   pts = pysesa_read.lasread(infile)

   Parameters
   ------------
   infile : str
   	LAS or LAZ file containing Nx3 point cloud

   Returns
   ----------
   data: ndarray
   	Nx3 point cloud, 32 bit precision

   '''

   inFile = File(infile, mode='r')

   return np.c_[ascol(scaled_x_dimension(inFile)), ascol(scaled_y_dimension(inFile)), ascol(scaled_z_dimension(inFile))].astype('float32')


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


   if len(data[0].decode().split(','))>1: # the file is comma delimited
      out = np.array([x.decode().split(',',2)[0:3] for x in data], dtype='float32')
      return out
   else:
      try: # space delimited
         out = np.array([x.decode().split(' ',2)[0:3] for x in data], dtype='float32')
         return out
      except: #tab delimited
         out = np.array([x.decode().split('\t',2)[0:3] for x in data], dtype='float32')
         return out
