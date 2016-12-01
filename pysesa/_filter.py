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
#    _____ ____
#   / __(_) / /____  _____
#  / /_/ / / __/ _ \/ ___/
# / __/ / / /_/  __/ /
#/_/ /_/_/\__/\___/_/

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

# import libraries
from __future__ import division
import numpy as np

try:
   from pykdtree.kdtree import KDTree
   pykdtree=1
except:
   #print "install pykdtree for faster kd-tree operations: https://github.com/storpipfugl/pykdtree"
   from scipy.spatial import cKDTree as KDTree
   pykdtree=0

# suppress divide and invalid warnings
import warnings
warnings.filterwarnings("ignore")
np.seterr(divide='ignore')
np.seterr(invalid='ignore')
np.seterr(over='ignore')
np.seterr(under='ignore')

##======================================================
def filt_stdev(coords, k=3, std_dev=2):

   kDTree = KDTree(coords, leafsize = 5)

   if pykdtree==1:
      dx, idx_knn = kDTree.query(coords[:, :], k = k)
   else:
      dx, idx_knn = kDTree.query(coords[:, :], k = k, n_jobs=-1)

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
