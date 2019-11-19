.. _pysesa.read:

pysesa.read module
===================

Custom fast (up to 3.5x faster than numpy's genfromtxt) txt file to numpy array
accepts comma, tab or space delimited files of 3 columns: x, y, and amplitude


Syntax
----------

You call the function like this::

  pts = pysesa_read.txtread(infile)

Parameters
------------
infile : str
	3-column ASCII file containing Nx3 point cloud

Returns
----------
data: ndarray
 	Nx3 point cloud, 32 bit precision

  .. image:: _static/pysesa_colour.jpg
