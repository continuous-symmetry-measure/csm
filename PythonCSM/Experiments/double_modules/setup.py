from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules= cythonize(
        [Extension(
            "module1",
            ["module1.pyx"],
            language='c++',
        ),
        Extension(
            "module2",
            ["module2.pyx"],
            language='c++',
        )])
)