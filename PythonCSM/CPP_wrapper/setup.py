__author__ = 'zmbq'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys

BOOST_ROOT = r'\boost\1_57_0\lib64-msvc-10.0' # Default, Windows only. Override in local_settings for now

try:
    from local_settings import *
except:
    pass

extra_compile_args = []
if sys.platform=='win32':
    library_dirs = ['../../openbabel-files/Windows/lib/x64/Release', '../../CSM/cmake/Release', BOOST_ROOT]
    libraries=['openbabel-2', 'csmlib']
elif sys.platform in ['linux', 'linux2']:
    library_dirs = ['../../openbabel-files/unix/lib', '../../CSM/cmake']
    libraries=['openbabel', 'csmlib', 'boost_log-mt', 'boost_log_setup-mt', 'boost_system-mt', 'boost_thread-mt', 'boost_filesystem-mt', 'boost_date_time-mt', 'pthread'],
    extra_compile_args = ['-fPIC']
elif sys.platform=='darwin':
    library_dirs = ['../../openbabel-files/unix/lib', '../../CSM/cmake']
    libraries=['openbabel', 'csmlib', 'boost_log', 'boost_log_setup', 'boost_system', 'boost_thread', 'boost_filesystem', 'boost_date_time', 'pthread'],

setup(
    ext_modules = cythonize(
        [Extension(
            "csm",
            ["csm.pyx"],
            language='c++',
            include_dirs=['../../CSM'],
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=extra_compile_args)]
    )
)
