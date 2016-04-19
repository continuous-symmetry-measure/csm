__author__ = 'zmbq'

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
import numpy
import os

FAST_CPPUTILS_DIR = "../../FastCPPUtils"
EIGEN_INCLUDE_DIR = "../../../include"

extra_compile_args = []
extra_link_args = []
if sys.platform == 'win32':
    extra_compile_args = ['/Ox']
    # extra_link_args = ['/debug']
elif sys.platform in ['linux', 'linux2']:
    extra_compile_args = ['-fPIC']

setup(
    name='fast',
    version='0.1',
    description='Fast C++ calculations',
    author='Chelem',
    author_email='contact@chelem.co.il',
    ext_modules=cythonize(
        [Extension(
            "fast",
            ["fast.pyx",
             os.path.join(FAST_CPPUTILS_DIR, "rpoly.c"),
             os.path.join(FAST_CPPUTILS_DIR, "math_wrappers.cpp")],
            language='c++',
            include_dirs=[numpy.get_include(), EIGEN_INCLUDE_DIR, FAST_CPPUTILS_DIR],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args)]
    )
)
