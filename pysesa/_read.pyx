"""
                        __
   ________  ____ _____/ /
  / ___/ _ \/ __ `/ __  / 
 / /  /  __/ /_/ / /_/ /  
/_/   \___/\__,_/\__,_/   
                          
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
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray txtread(infile):
   '''
   custom fast (up to 3.5x faster than numpy's genfromtxt) txt file to numpy array
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
         
