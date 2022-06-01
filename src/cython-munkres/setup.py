from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
import numpy

setup(
    name='cython_munkres',
    url='https://github.com/jfrelinger/cython-munkres-wrapper',
    ext_modules = cythonize(
        [Extension("cython_munkres",
                   ["src/cython_munkres.pyx", "src/cpp/Munkres.cpp"],
                   include_dirs = [numpy.get_include(), 'src/cpp'],
                   language='c++')]),
    version = '1.0',
    description='Munkres implemented in c++ wrapped by cython',
    author='Jacob Frelinger',
    author_email='jacob.frelinger@duke.edu',
    requires=['numpy (>=1.3.0)', 'cython (>=0.15.1)']
)

