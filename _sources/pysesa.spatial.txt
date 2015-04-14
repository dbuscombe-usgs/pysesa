.. _pysesa.spatial:

pysesa.spatial module
======================

Calculate spatial statistics of a Nx3 point cloud

Syntax
----------

You call the function like this::

  stats = pysesa.spatial(points).getdata()

or::

  centroids = pysesa.spatial(points).getcentroid()

or::

  stats = pysesa.spatial(points).getstats()

Parameters
-------------
points : ndarray
	Nx3 point cloud

Returns [requested through .getdata()]
---------------------------------------
self.data: list
	x = centroid in horizontal coordinate

        y = centroid in laterial coordinate

        z_mean = centroid in amplitude

        z_max = max amplitude

        z_min = min amplitude

        z_range = range in amplitude

        sigma = standard deviation of amplitudes

        skewness = skewness of amplitudes

        kurtosis = skewness of amplitudes

        n = number of 3D coordinates


Returns [requested through .getcentroid()]
--------------------------------------------
self.centroid: list
   	1x3 point cloud centroid [x,y,z]

Returns [requested through .getstats()]
----------------------------------------
self.stats: list
        z_mean = centroid in amplitude

        z_max = max amplitude

        z_min = min amplitude

        z_range = range in amplitude

        sigma = standard deviation of amplitudes

        skewness = skewness of amplitudes

        kurtosis = skewness of amplitudes

        n = number of 3D coordinates

  .. image:: _static/pysesa_colour.jpg
