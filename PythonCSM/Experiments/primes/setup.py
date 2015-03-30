from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='cython_primes',
    ext_modules = cythonize("primes.pyx")
)