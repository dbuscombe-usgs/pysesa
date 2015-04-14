.. _pysesa.partition:

pysesa.partition module
========================

Partition a Nx3 point cloud into M windows of nx3 points
with specified spacing between centroids of adjacent windows
and with specified overlap between windows.
Implemented using a binary search tree for fast nearest neighbour 
point check with boundary pruning

Syntax
----------

You call the function like this::

  nr_pts = pysesa.partition(toproc, out, res, mxpts, minpts, prc_overlap).getdata()

Parameters
-----------
toproc : ndarray
	Nx3 point cloud

Other Parameters
-----------------
out : float, *optional* [default = 0.5]
	output grid resolution
res : float, *optional* [default = 0.05]
	spatial grid resolution to create a grid for the boundary pruning
mxpts : float, *optional* [default = 1024]
	maximum number of points allowed in a window
minpts : float, *optional* [default = 16]
	minimum number of points allowed in a window
prc_overlap : float, *optional"  [default = 0]
     percentage overlap between windows

Returns
----------
self.data: list
	list of M ndarrays, each containing n indices 
     of original point cloud, toproc, to partition space 
     to create M windows

  .. image:: _static/pysesa_colour.jpg
