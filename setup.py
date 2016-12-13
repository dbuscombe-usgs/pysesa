#!/usr/bin/env python
# -*- coding: utf-8 -*-
## PySESA (Python program for Spatially Explicit Spectral Analysis) 
## has been developed at the Grand Canyon Monitorinf & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/pysesa>
##
##This software is in the public domain because it contains materials that originally came from 
##the United States Geological Survey, an agency of the United States Department of Interior. 
##For more information, see the official USGS copyright policy at 
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
PySESA - a Python framework for Spatially Explicit Spectral Analysis

PySESA is an open-source project dedicated to provide a generic Python framework 
for spatially explicit statistical analyses of point clouds and other geospatial data, 
in the spatial and frequency domains, for use in the geosciences

The program is detailed in:
   Buscombe, D. (2016) "Computational considerations for spatially explicit spectral analysis of point clouds and geospatial data", 86, 92-108, 10.1016/j.cageo.2015.10.004.

:Author:  
    Daniel Buscombe
    Grand Canyon Monitoring and Research Center
    United States Geological Survey
    Flagstaff, AZ 86001
    dbuscombe@usgs.gov

For more information visit http://dbuscombe-usgs.github.io/pysesa/

:install:
    pip install pysesa
    python setup.py install
    sudo python setup.py install
    
:test:
    python -c "import pysesa; pysesa.test()"

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
           "Please install numpy first, it is needed before installing PySESA.")
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
      from Cython.Build import cythonize
   except:
      msg = ("No module named Cython. "
           "Please install Cython first, it is needed before installing pysesa.")
      raise ImportError(msg)

# Read version from distmesh/__init__.py
with open(os.path.join('pysesa', '__init__.py')) as f:
    line = f.readline()
    while not line.startswith('__version__'):
        line = f.readline()
exec(line, globals())

ext_modules = [ ]
cmdclass = { }

#ext_modules += cythonize(Extension('_loadtxt',sources=['pysesa/loadtxt.pyx', 'pysesa/loadtxt.cpp'], language="c++"))

if USE_CYTHON:
    ext_modules += [
        Extension("pysesa.plot", [ "pysesa/plot.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.read", [ "pysesa/read.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.write", [ "pysesa/_write.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.detrend", [ "pysesa/detrend.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.lengthscale", [ "pysesa/lengthscale.pyx" ],
        include_dirs=[np.get_include()]),        
        Extension("pysesa.partition", [ "pysesa/partition.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spectral", [ "pysesa/spectral.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spatial", [ "pysesa/_spatial.pyx" ],
        include_dirs=[np.get_include()]),    
        Extension("pysesa.sgolay", [ "pysesa/_sgolay.pyx" ],
        include_dirs=[np.get_include()]),
        Extension('_RunningStats',sources=['pysesa/RunningStats_wrap.cxx', 'pysesa/RunningStats.cpp']),
        ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("pysesa.plot", [ "pysesa/plot.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.read", [ "pysesa/read.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.write", [ "pysesa/_write.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.detrend", [ "pysesa/detrend.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.lengthscale", [ "pysesa/lengthscale.c" ],
        include_dirs=[np.get_include()]), 
        Extension("pysesa.partition", [ "pysesa/partition.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spectral", [ "pysesa/spectral.c" ],
        include_dirs=[np.get_include()]),
        Extension("pysesa.spatial", [ "pysesa/_spatial.c" ],
        include_dirs=[np.get_include()]),    
        Extension("pysesa.sgolay", [ "pysesa/_sgolay.c" ],
        include_dirs=[np.get_include()]),
        Extension('_RunningStats',sources=['pysesa/RunningStats_wrap.cxx', 'pysesa/RunningStats.cpp']),
    ]
install_requires = [
    'numpy','scipy','matplotlib', 'cython', 'statsmodels', 'ift_nifty', 'joblib','dask', 'toolz', 'dill', 'laspy'
]

def setupPackage():
   setup(name='pysesa',
         version=__version__,
         description='PySESA is an open-source project dedicated to provide a generic Python framework for spatially explicit statistical analyses of point clouds and other geospatial data, in the spatial and frequency domains, for use in the geosciences. The program is detailed in  Buscombe, D. (2016) "Computational considerations for spatially explicit spectral analysis of point clouds and geospatial data", 86, 92-108, 10.1016/j.cageo.2015.10.004.',
         #long_description=long_description,
         classifiers=[
             'Intended Audience :: Science/Research',
             'Intended Audience :: Developers',
             'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
             'Programming Language :: Python',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Cython',
             'Topic :: Scientific/Engineering',
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
         package_data={'pysesa': ['*.xyz','*.h']}
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

