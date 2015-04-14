.. _pysesa.lengthscale:

pysesa.lengthscale module
==========================

Calculates the integral lengthscale of a Nx3 point cloud
using 1 of 3 available methods
and also returns the tapered 2D grid of 3D pointcloud for spectral analysis

Syntax
----------

You call the function like this::

   im = pysesa.lengthscale(points, res, lentype, taper, method).getdata()

or::

   lengthscale = pysesa.lengthscale(points, res, lentype, taper, method).getlengthscale()

Parameters
------------
points : ndarray
	Nx3 point cloud

Other Parameters
------------------
res : float, *optional* [default = 0.05]
     spatial grid resolution to create a grid
lentype : int, *optional* [default = 0, l<0.5]
	lengthscale type:
     1, l<0.5

     2, l<1/e

     3, l<0

taper : int, *optional* [default = Hanning]
	flag for taper type:
     1, Hanning (Hann)

     2, Hamming

     3, Blackman

     4, Bartlett

method : str, *optional* [default = 'nearest']
	gridding type

Returns [requested through .getdata()]
-----------------------------------------
self.data: ndarray
	tapered 2D grid of 3D pointcloud

Returns [requested through .getlengthscale()]
----------------------------------------------
self.lengthscale: float
	integral lengthscale

  .. image:: _static/pysesa_colour.jpg


