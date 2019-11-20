# 8""""8        8""""8 8"""" 8""""8 8""""8
# 8    8 e    e 8      8     8      8    8
# 8eeee8 8    8 8eeeee 8eeee 8eeeee 8eeee8
# 88     8eeee8     88 88        88 88   8
# 88       88   e   88 88    e   88 88   8
# 88       88   8eee88 88eee 8eee88 88   8

## python setup.py build_ext --inplace
from distutils.core import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy as np
from Cython.Build import cythonize


extensions = [Extension("read", ["src/_read.pyx"], include_dirs=[np.get_include()])]
setup(
    name="read",
    ext_modules=cythonize(extensions),
)

extensions = [Extension("detrend", ["src/_detrend.pyx"], include_dirs=[np.get_include()])]
setup(
    name="detrend",
    ext_modules=cythonize(extensions),
)

extensions = [Extension("lengthscale", ["src/_lengthscale.pyx"], include_dirs=[np.get_include()])]
setup(
    name="lengthscale",
    ext_modules=cythonize(extensions),
)


extensions = [Extension("spatial", ["src/_spatial.pyx"], include_dirs=[np.get_include()])]
setup(
    name="spatial",
    ext_modules=cythonize(extensions),
)


extensions = [Extension("spectral", ["src/_spectral.pyx"], include_dirs=[np.get_include()])]
setup(
    name="spectral",
    ext_modules=cythonize(extensions),
)


extensions = [Extension("write", ["src/_write.pyx"], include_dirs=[np.get_include()])]
setup(
    name="write",
    ext_modules=cythonize(extensions),
)


# extensions = [Extension("plot", ["src/_plot.pyx"], include_dirs=[np.get_include()])]
# setup(
#     name="plot",
#     ext_modules=cythonize(extensions),
# )


extensions = [Extension("partition", ["src/_partition.pyx"], include_dirs=[np.get_include()])]
setup(
    name="partition",
    ext_modules=cythonize(extensions),
)
