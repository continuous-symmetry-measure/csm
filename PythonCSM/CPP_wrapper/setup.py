__author__ = 'zmbq'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(
        [Extension(
            "csm",
            ["csm.pyx"],
            language='c++',
            include_dirs=['../../CSM'],
            library_dirs=['../../openbabel-files/Windows/lib/x64/Release', '../../CSM/cmake/Release', r'd:\boost\1_57_0\lib64-msvc-10.0'],
            libraries=['openbabel-2', 'csmlib'],
        )])
)