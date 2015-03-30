# encoding: utf-8
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

__version__ = '0.0.1'

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from pysesa._main import pysesa
from pysesa.test import *


