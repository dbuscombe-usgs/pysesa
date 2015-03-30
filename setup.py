#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pysesa - a Python framework for .......

pysesa is an open-source project dedicated to provide a Python framework for ...

For more information visit http://dbuscombe-usgs.github.io/pysesa/

:install:
    python setup.py install
    sudo python setup.py install
    
:test:
    python -c "import pysesa; pysesa.test.dotest()"

:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
    
    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information, 
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only 
    and does not imply endorsement by the U.S. government.
    
"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys, glob
import inspect
try:
    import numpy as np
except:
    msg = ("No module named numpy. "
           "Please install numpy first, it is needed before installing pysesa.")
    raise ImportError(msg)

from distutils.core import setup
from distutils.extension import Extension
#from setuptools import setup, Extension

# Directory of the current file 
SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))

# Set this to True to enable building extensions using Cython.
# Set it to False to build extensions from the C file (that
# was previously created using Cython).
# Set it to 'auto' to build with Cython if available, otherwise
# from the C file.
USE_CYTHON = True

if USE_CYTHON:
   try:
      from Cython.Distutils import build_ext
   except:
      msg = ("No module named Cython. "
           "Please install Cython first, it is needed before installing PyHum.")
      raise ImportError(msg)

#if USE_CYTHON:
#    try:
#        from Cython.Distutils import build_ext
#    except ImportError:
#        if USE_CYTHON=='auto':
#            USE_CYTHON=False
#        else:
#            raise

# Read version from distmesh/__init__.py
with open(os.path.join('pysesa', '__init__.py')) as f:
    line = f.readline()
    while not line.startswith('__version__'):
        line = f.readline()
exec(line, globals())

ext_modules = [ ]
cmdclass = { }

if USE_CYTHON:
    ext_modules += [
        Extension("pysesa.partition", [ "pysesa/_partition.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spec", [ "pysesa/_spec.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.worker", [ "pysesa/_worker.pyx" ],
        include_dirs=[np.get_include()]),    
        Extension("pysesa.sgolay", [ "pysesa/_sgolay.pyx" ],
        include_dirs=[np.get_include()]),
        Extension('_RunningStats',sources=['pysesa/RunningStats_wrap.cxx', 'pysesa/RunningStats.cpp']),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("pysesa.partition", [ "pysesa/_partition.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spec", [ "pysesa/_spec.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.worker", [ "pysesa/_worker.c" ],
        include_dirs=[np.get_include()]),    
        Extension("pysesa.sgolay", [ "pysesa/_sgolay.c" ],
        include_dirs=[np.get_include()]),
        Extension('_RunningStats',sources=['pysesa/RunningStats_wrap.cxx', 'pysesa/RunningStats.cpp']),
    ]
install_requires = [
    'numpy','scipy','matplotlib', 'cython', 
]
#long_description = open('README.md').read()

def setupPackage():
   setup(name='pysesa',
         version=__version__,
         description='Python/Cython scripts to ..... detailed in Buscombe, "Computational considerations for spatially explicit spectral analysis of point clouds", forthcoming.',
         #long_description=long_description,
         classifiers=[
             'Intended Audience :: Science/Research',
             'Intended Audience :: Developers',
             'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
             'Programming Language :: Python',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Cython',
             'Topic :: Scientific/Engineering',
             'Topic :: Scientific/Engineering :: Physics',
         ],
         keywords='point clouds',
         author='Daniel Buscombe',
         author_email='dbuscombe@usgs.gov',
         url='https://github.com/dbuscombe-usgs/pysesa',
         download_url ='https://github.com/dbuscombe-usgs/pysesa/archive/master.zip',
         install_requires=install_requires,
         license = "GNU GENERAL PUBLIC LICENSE v3",
         packages=['pysesa'],
         cmdclass = cmdclass,
         ext_modules=ext_modules,
         platforms='OS Independent',
         package_data={'pysesa': ['*.xyz']}
   )

if __name__ == '__main__':
    # clean --all does not remove extensions automatically
    if 'clean' in sys.argv and '--all' in sys.argv:
        import shutil
        # delete complete build directory
        path = os.path.join(SETUP_DIRECTORY, 'build')
        try:
            shutil.rmtree(path)
        except:
            pass
        # delete all shared libs from lib directory
        path = os.path.join(SETUP_DIRECTORY, 'pysesa')
        for filename in glob.glob(path + os.sep + '*.pyd'):
            try:
                os.remove(filename)
            except:
                pass
        for filename in glob.glob(path + os.sep + '*.so'):
            try:
                os.remove(filename)
            except:
                pass
    setupPackage()

