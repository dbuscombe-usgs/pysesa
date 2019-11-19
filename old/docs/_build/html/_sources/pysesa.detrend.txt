.. _pysesa.pysesa:

pysesa.detrend module
======================

Detrend a Nx3 point cloud by specified method

Syntax
----------

You call the function like this::

  detrended_pts = pysesa.detrend(points, proctype, res, method).getdata()

Parameters
----------
points : ndarray
	Nx3 point cloud
proctype : int
	type of detrending.
     1 = remove mean

     2 = remove Ordinary least squares plane

     3 = remove Robust linear model plane

     4 = remove Orthogonal Distance Regression plane

     5 = remove Savitsky-Golay digital filter, order 1

Other Parameters
----------
res : float, *optional* [default = 0.05]
	for proctype==4 only
     spatial grid resolution to create a grid
method : str, *optional* [default = 'nearest']
	for proctype==4 only
	gridding type

Returns
----------
self.data: ndarray
	Nx3 detrended point cloud

  .. image:: _static/pysesa_colour.jpg
