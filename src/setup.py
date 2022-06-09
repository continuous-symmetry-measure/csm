from setuptools import setup #if this import statement dooes not come first, weird errors happen
#see: https://stackoverflow.com/questions/21594925/error-each-element-of-ext-modules-option-must-be-an-extension-instance-or-2-t
from distutils.extension import Extension
import setuptools
import sys
import numpy
import os
import re
import glob
from setuptools.command.build_ext import build_ext as _build_ext
import install_requirements

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
EIGEN_INCLUDE_DIR_1 = "include"
EIGEN_INCLUDE_DIR_2 = "include/Eigen"

CPP_WRAPPER_DIR = "csm/CPP_wrapper"

extra_compile_args = []
extra_link_args = []
if sys.platform == 'win32':
    extra_compile_args = ['/Ox']
    # extra_link_args = ['/debug']
elif sys.platform in ['linux', 'linux2']:
    extra_compile_args = ['-fPIC']

def get_version():
    pathname = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csm', 'version.py')
    with open(pathname, "r") as ver_inp:
        code = ver_inp.read()
    pattern = r"__version__\s*=\s*'(?P<version>.*)'"
    match = re.match(pattern, code)
    if not match:
        raise ValueError("Version file must contain one line: __version__='...'")
    version = match.group("version")
    return version


class PrepareCommand(setuptools.Command):
    description = "Build fast.pyx so there's no cython dependence in installation"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("running prepare command")
        self.copy_source_files()
        self.convert_to_c()

    def copy_source_files(self):
        #this may be used for copying openbabel file eventually?
        pass

    def convert_to_c(self):
        #creates fast.h and fast.c in cpp_wrapper folder
        print('Converting csm pyx files to C++ sources...')
        pyx = './csm/CPP_wrapper/fast.pyx'
        self.cython(pyx)

    def cython(self, pyx):
        from Cython.Compiler.CmdLine import parse_command_line
        from Cython.Compiler.Main import compile
        options, sources = parse_command_line(['-2', '-v', '--cplus', pyx])
        result = compile(sources, options)
        if result.num_errors > 0:
            print('Errors converting %s to C++' % pyx, file=sys.stderr)
            raise Exception('Errors converting %s to C++' % pyx)
        self.announce('Converted %s to C++' % pyx)

csm_version = get_version()
openbabel_dependency = install_requirements.get_openbabel_dependency()
print("Packaging CSM version %s" % csm_version)
setup(
    name='csm',
    version=csm_version,
    packages=['csm.calculations', 'csm.calculations.approx', 'csm.input_output', 'csm.molecule', 'csm.main', 'csm',],
    setup_requires=['numpy>=1.10'],
    install_requires=['numpy>=1.10', openbabel_dependency, 'scipy>=1.7.3'],
    include_package_data=True,
    license='Chelem',  # example license
    description='The Continuous Symmetry Measure',
    long_description=README,
    url='http://www.csm.huji.ac.il/new/',
    author='The Research Software Company',
    author_email='itay@chelem.co.il',

    # The csm command
    entry_points={
        'console_scripts': [
            'csm = csm.main.csm_run:run_no_return',
            'norm_csm = csm.main.normcsm:run_norm_no_return',
            'csmsymm = csm.main.csmsymm:run_no_return'
        ]
    },

    classifiers=[
        'Development Status:: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],

    cmdclass={
        'prepare': PrepareCommand,
    },

    # The Cython extension module
    # We do not use Cython itself but the Cython output files
    # Note that header files that are required by the compilation are specified in MANIFEST.in
    ext_modules=[
        Extension(
            "csm.fast",
            [os.path.join(CPP_WRAPPER_DIR, "fast.cpp"),
             os.path.join(FAST_CPPUTILS_DIR, "rpoly.c"),
             os.path.join(FAST_CPPUTILS_DIR, "math_wrappers.cpp"),
             ],
            language='c++',
            include_dirs=[EIGEN_INCLUDE_DIR_1, EIGEN_INCLUDE_DIR_2, FAST_CPPUTILS_DIR, numpy.get_include()],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args),
    ],


)
