__author__ = 'zmbq'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(
        [Extension(
            "cpluspluslib",
            ["cpluspluslib.pyx"],
            language='c++',
            include_dirs=['../cpluspluslib'],
            library_dirs=['../cpluspluslib/x64/debug'],
            libraries=['cpluspluslib'],
            extra_compile_args=["/Zi", "/Od"],
            extra_link_args=["/debug"],
        )])
)