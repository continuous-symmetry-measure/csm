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
            library_dirs=['../../openbabel-files/Windows/lib/x64/Release', '../../CSM/cmake/Release', r'c:\boost\1_57_0\lib64-msvc-10.0'],
            libraries=['openbabel-2', 'csmlib'],
#            libraries=['openbabel', 'csmlib', 'boost_log-mt', 'boost_log_setup-mt', 'boost_system-mt', 'boost_thread-mt', 'boost_filesystem-mt', 'boost_date_time-mt', 'pthread'],
#            libraries=['openbabel', 'csmlib', 'boost_log', 'boost_log_setup', 'boost_system', 'boost_thread', 'boost_filesystem', 'boost_date_time', 'pthread'],
#            library_dirs=['../../openbabel-files/unix/lib', '../../CSM/cmake'],
#	    extra_compile_args=['-fPIC'],
        )])
)
