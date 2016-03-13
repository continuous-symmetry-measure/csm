__author__ = 'zmbq'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
import numpy

try:
    from local_settings import *
except:
    pass

extra_compile_args = []
extra_link_args = []
if sys.platform == 'win32':
    extra_compile_args = ['/Ox']  # Debug info, no optimization
elif sys.platform in ['linux', 'linux2']:
    extra_compile_args = ['-fPIC', '-O3']
elif sys.platform == 'darwin':
    extra_compile_args = ['-O3']

setup(
    ext_modules=cythonize(
        [Extension(
            "*",
            ["playground.pyx"],
            language='c++',
            include_dirs=[numpy.get_include()],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args)]
    )
)
