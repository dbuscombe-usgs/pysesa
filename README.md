### Pysesa

Python/Cython scripts to: 

1. one
2. two
3. three

<!--![alt tag](http://dbuscombe-usgs.github.io/figs/class_R01560.png)-->
<!--*Sand dunes on the bed of the Colorado River in Grand Canyon*-->

### Contributing & Credits

 Author |    Daniel Buscombe 
 ------ | ---------------
        |  Grand Canyon Monitoring and Research Center
        | United States Geological Survey
        | Flagstaff, AZ 86001
        | dbuscombe@usgs.gov
Version: 0.0.1    |  Revision: Apr, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of pysesa software
This software is in the public domain because it contains materials that originally came 
from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 

```
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
```

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government. 

### Contents

The programs in this package are as follows:

1. py_sesa

2. py_sesa_detrend

3. py_sesa_partition

4. py_sesa_spectral

5. py_sesa_spatial

6. py_sesa_plot

These are all command-line programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options


## Setup

### Automatic Installation from PyPI 

```
pip uninstall pysesa (removes previous installation)
pip install pysesa
```

Automatic Installation from github:

```
git clone git@github.com:dbuscombe-usgs/pysesa.git
cd pysesa
python setup.py install
```

or a local installation:

```
python setup.py install --user
```

or with admin privileges, e.g.:

```
sudo python setup.py install
```

### Notes for Windows/Anaconda users

Assuming a Anaconda distribution which comes with almost all required program dependencies:

```
pip install ??????
pip uninstall pysesa (removes previous installation)
pip install pysesa
```


### Notes for Linux users

You could try before you install, using a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install numpy
pip install cython
pip install scipy
pip install joblib
pip install statsmodels
pip install matplotlib
pip install pysesa
python -c "import pysesa; pysesa.test.dotest()"
deactivate (or source venv/bin/deactivate)
```

The results will live in "venv/lib/python2.7/site-packages/pysesa"


###Manual Installation

PYTHON LIBRARIES YOU MAY NEED TO INSTALL TO USE pysesa:

1. Nifty ()
2. SciPy (http://www.scipy.org/scipylib/download.html)
3. Numpy (http://www.scipy.org/scipylib/download.html)
4. Matplotlib (http://matplotlib.org/downloads.html)
5. cython ()
6. statsmodels ()

All of the above are available through pip (https://pypi.python.org/pypi/pip) and easy_install (https://pythonhosted.org/setuptools/easy_install.html)

### Test

A test can be carried out by running the supplied script:

```
python -c "import pysesa; pysesa.test.dotest()"
```

which carries out the following operations:

```


```

on the following files:

example_100000.xyz

and results in 

### Support

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 



