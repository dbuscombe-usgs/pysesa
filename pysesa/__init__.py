# encoding: utf-8
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

__version__ = '0.0.31'

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import matplotlib
matplotlib.use('TKagg')

from pysesa._pysesa import process
from pysesa.test import *

import pysesa.read
import pysesa.write
from pysesa.partition import partition
from pysesa.sgolay import sgolay
from pysesa.spatial import spatial
from pysesa.spectral import spectral
from pysesa.lengthscale import lengthscale
from pysesa.detrend import detrend

from pysesa.plot import plot

import pysesa._filter as filter
