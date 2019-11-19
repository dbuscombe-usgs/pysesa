.. _pysesa.sgolay:

pysesa.sgolay module
=====================

Create a Savitsky-Golay digital filter from a 2D signal
based on code from http://www.scipy.org/Cookbook/SavitzkyGolay

Syntax
----------

You call the function like this::

  Z = pysesa.sgolay(z, window_size, order).getdata()

Parameters
----------
z : array_like, shape (N,)
  the 2D signal.
window_size : int
   the length of the window. Must be an odd integer number.
order : int
   the order of the polynomial used in the filtering.
   Must be less than `window_size` - 1.

Returns
-------
self.data : ndarray, shape (N)
   the smoothed signal.

References
----------
.. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
   Data by Simplified Least Squares Procedures. Analytical
   Chemistry, 1964, 36 (8), pp 1627-1639.
.. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
   W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
   Cambridge University Press ISBN-13: 9780521880688

  .. image:: _static/pysesa_colour.jpg
