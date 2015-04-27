.. _pysesa.test:

pysesa.test module
==================

You call the function like this::
  python -c "import pysesa; pysesa.test.dotest()"

which calls pysesa.process with the following default parameters::

  # copy files over to somewhere read/writeable
  dircopy(pysesa.__path__[0], os.path.expanduser("~")+os.sep+'pysesa_test')
  shutil.copy(pysesa.__path__[0]+os.sep+'example_100000pts.xyz', os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'test.DAT')

  # general settings   
  infile = os.path.expanduser("~")+os.sep+'pysesa_test'+os.sep+'example_100000pts.xyz' 

  out = 0.5 #m output grid
  detrend = 4 #ODR plane
  proctype = 1 #Processing spectral parameters (no smoothing)
  mxpts = 512 # max pts per window
  res = 0.05 #cm internal grid resolution
  nbin = 20 #number of bins for spectral binning
  lentype = 1 # l<0.5
  taper = 1 # Hann taper
  prc_overlap = 0 # no overlap between successive windows
  minpts = 16 # min pts per window

  pysesa.process(infile, out, detrend, proctype, mxpts, res, nbin, lentype, minpts, taper, prc_overlap)

  .. image:: _static/pysesa_colour.jpg
