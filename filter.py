# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8

# import libraries
from __future__ import division
import numpy as np
from sklearn.neighbors import KDTree

# suppress divide and invalid warnings
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

##======================================================
def filt_stdev(coords, k=3, std_dev=2):

   kDTree = KDTree(coords, leaf_size = 5) #leafsize = 5)
   dx, idx_knn = kDTree.query(coords[:, :], k = k)

   dx, idx_knn = dx[:,1:], idx_knn[:,1:]

   distances = np.sum(dx, axis=1)/(k - 1.0)
   valid_distances = np.shape(distances)[0]

   #Estimate the mean and the standard deviation of the distance vector
   sum = np.sum(distances)
   sq_sum = np.sum(distances**2)

   mean = sum / float(valid_distances)
   variance = (sq_sum - sum * sum / float(valid_distances)) / (float(valid_distances) - 1)
   stddev = np.sqrt (variance)

   # a distance that is bigger than this signals an outlier
   distance_threshold = mean + std_dev * stddev
   idx = np.nonzero(distances < distance_threshold)

   return idx, np.copy(coords[idx])


   # if pykdtree==1:
   #    dx, idx_knn = kDTree.query(coords[:, :], k = k)
   # else:
   #    dx, idx_knn = kDTree.query(coords[:, :], k = k, n_jobs=-1)


#from pykdtree.kdtree import KDTree
#pykdtree=1
# except:
#    #print "install pykdtree for faster kd-tree operations: https://github.com/storpipfugl/pykdtree"
#    from scipy.spatial import cKDTree as KDTree
#    pykdtree=0
