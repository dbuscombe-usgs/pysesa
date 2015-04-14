.. _pysesa.write:

pysesa.write module
====================

Custom fast numpy array to comma-delimited ASCII txt file

Syntax
----------

You call the function like this::

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

  .. image:: _static/pysesa_colour.jpg
