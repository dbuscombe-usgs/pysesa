.. _getting_started:


***************
Getting started
***************

.. _about:

About
======

PySESA - a Python framework for Spatially Explicit Spectral Analysis

PySESA is an open-source project dedicated to provide a generic Python framework 
for spatially explicit statistical analyses of point clouds and other geospatial data, 
in the spatial and frequency domains, for use in the geosciences

The program is detailed in:
Buscombe, D. "Spatially explicit spectral analysis of point clouds and geospatial data", forthcoming.

For the source code visit `my project github site <http://dbuscombe-usgs.github.io/pysesa/>`_


.. _license:

License
========

This software is in the public domain because it contains materials that
originally came from the United States Geological Survey, an agency of the
United States Department of Interior. For more information, 
see `the official USGS copyright policy <http://www.usgs.gov/visual-id/credit_usgs.html#copyright>`_

Any use of trade, product, or firm names is for descriptive purposes only 
and does not imply endorsement by the U.S. government.

This software is issued under the `GNU Lesser General Public License, Version 3 <http://www.gnu.org/copyleft/lesser.html>`_


.. _setup:

Setup
========

Automatic Installation from PyPI::


  pip uninstall pysesa (removes previous installation)
  pip install pysesa


Automatic Installation from github::


  git clone git@github.com:dbuscombe-usgs/pysesa.git
  cd pysesa
  python setup.py install


or a local installation::


  python setup.py install --user


or with admin privileges, e.g.::


  sudo python setup.py install


.. _virtualenv:

Virtual environment
====================

You could try before you install, using a virtual environment::

  virtualenv venv
  source venv/bin/activate
  pip install numpy
  pip install cython
  pip install scipy
  pip install joblib
  pip install statsmodels
  pip install matplotlib
  pip install ift_nifty
  pip install pysesa
  python -c "import pysesa; pysesa.test()"
  deactivate #(or source venv/bin/deactivate)

The results will live in "venv/lib/python2.7/site-packages/pysesa"


.. _manualinstall:

Manual installation
====================

Python libraries you need to have installed to use pysesa:

1. `Nifty <http://www.mpa-garching.mpg.de/ift/nifty/index.html>`_
2. `SciPy <http://www.scipy.org/scipylib/download.html>`_
3. `Numpy <http://www.scipy.org/scipylib/download.html>`_
4. `Matplotlib <http://matplotlib.org/downloads.html>`_
5. `cython <http://cython.org/>`_
6. `statsmodels <http://statsmodels.sourceforge.net/>`_

All of the above are available through `pip <https://pypi.python.org/pypi/pip>`_ and `easy_install <https://pythonhosted.org/setuptools/easy_install.html>`_


.. _test:

Test
======

A test can be carried out by running the supplied script::

  python -c "import pysesa; pysesa.test()"

which carries out the following operations::

  # general settings   
  infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_100000pts.xyz' 

  out = 1 #m output grid
  detrend = 4 #ODR plane
  proctype = 1 #Processing spectral parameters (no smoothing)
  mxpts = 1024 # max pts per window
  res = 0.05 #cm internal grid resolution
  nbin = 20 #number of bins for spectral binning
  lentype = 1 #l less than 0.5
  taper = 1 #Hann taper
  prc_overlap = 0 #no overlap between successive windows
  minpts = 64 #min pts per window

  pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)


.. _support:

Support
=========

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 

  .. image:: _static/pysesa_colour.jpg

