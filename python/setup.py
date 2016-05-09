import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import sys
import numpy

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))


# Fix the numpy bootstrap problem (http://stackoverflow.com/a/21621689/871910)
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())
        print("Include dirs are: ", self.include_dirs)

# Cython definitions
FAST_CPPUTILS_DIR = "FastCPPUtils"
EIGEN_INCLUDE_DIR = "../include"
CPP_WRAPPER_DIR = "csm/CPP_wrapper"

extra_compile_args = []
extra_link_args = []
if sys.platform == 'win32':
    extra_compile_args = ['/Ox']
    # extra_link_args = ['/debug']
elif sys.platform in ['linux', 'linux2']:
    extra_compile_args = ['-fPIC']


setup(
    name='csm',
    version='0.7.0',
    packages=['csm.calculations', 'csm.input_output', 'csm.molecule', 'csm.main'],
    setup_requires=['numpy>=1.10'],
    install_requires=['numpy>=1.10', 'openbabel>=1.8'],
    include_package_data=True,
    license='Chelem',  # example license
    description='The Continuous Symmetry Measure',
    long_description=README,
    url='http://www.csm.huji.ac.il/new/',
    author='The Research Software Company',
    author_email='itay@chelem.co.il',

    # Files required by the extension
    package_data={
        "csm.main": [os.path.join(FAST_CPPUTILS_DIR, "*.*"), 'README.md'],
    },

    # The Cython extension module
    ext_modules=[Extension(
        "csm.fast",
        [os.path.join(CPP_WRAPPER_DIR, "fast.cpp"),
         os.path.join(FAST_CPPUTILS_DIR, "rpoly.c"),
         os.path.join(FAST_CPPUTILS_DIR, "math_wrappers.cpp"),
         ],
        language='c++',
        include_dirs=[EIGEN_INCLUDE_DIR, FAST_CPPUTILS_DIR, numpy.get_include()],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args),],

    # The csm command
    entry_points={
        'console_scripts': [
            'csm = csm.main.csm_run:run_no_return',
        ]
    },

    classifiers=[
        'Development Status:: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
)
