__author__ = 'zmbq'

from distutils.core import setup
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
    library_dirs = ['../../openbabel-files/Windows/lib/x64/Release', '../../CSM/cmake/RelWithDebInfo', BOOST_ROOT]
    libraries = ['csmlib']
    # extra_compiler_args = ['/Zi', '/Od']  # Debug info, no optimization
    # extra_link_args = ['/debug']
elif sys.platform in ['linux', 'linux2']:
    library_dirs = ['../../openbabel-files/unix/lib', '../../CSM/cmake']
    libraries = ['csmlib', 'boost_log', 'boost_log_setup', 'boost_system', 'boost_thread', 'boost_filesystem', 'boost_date_time', 'pthread']
    extra_compile_args = ['-fPIC']
elif sys.platform == 'darwin':
    library_dirs = ['../../openbabel-files/unix/lib', '../../CSM/cmake']
    libraries = ['csmlib', 'boost_log-mt', 'boost_log_setup-mt', 'boost_system-mt', 'boost_thread-mt', 'boost_filesystem-mt', 'boost_date_time-mt', 'pthread']

setup(
    ext_modules=cythonize(
        [Extension(
            "*",
            ["permuters.pyx", "cache.pyx"],
            language='c++',
            include_dirs=[numpy.get_include()],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args)]
    )
)
