__author__ = 'zmbq'

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
import numpy

BOOST_ROOT = r'd:\boost\1_59_0\lib64' # Default, Windows only. Override in local_settings for now

try:
    from local_settings import *
except:
    pass

extra_compile_args = []
extra_link_args = []
if sys.platform == 'win32':
    extra_compile_args = ['/Ox']
    # extra_link_args = ['/debug']
elif sys.platform in ['linux', 'linux2']:
    extra_compile_args = ['-fPIC']

setup(
    name='CPP_Wrapper',
    version='0.1',
    description='Fast C++ calculations',
    author='Chelem',
    author_email='contact@chelem.co.il',
    ext_modules=cythonize(
        [Extension(
            "*",
            ["fast.pyx", "../../FastCPPUtils/rpoly.c", "../../FastCPPUtils/math_wrappers.cpp"],
            language='c++',
            include_dirs=[numpy.get_include(), '../../include', '../../FastCPPUtils'],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args)]
    )
)
