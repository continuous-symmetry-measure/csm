# Windows setup of OpenBabel
#
# Starting with Python 3.8, Python does not search for DLLs on the Path,
# so the openbabel extension can't locate the openbabel dlls. We need to explicitly
# add tell Python to look for the DLL in the right place, using os.add_dll_directory.
#
# This code needs to run before any openbabel imports
import sys
import os
from csm.input_output import formatters
from csm.input_output.formatters import csm_log as print

openbabel_prepared = False
def prepare_openbabel():
    formatters.csm_out_pipe = sys.stderr
    print("prepare_openbabel here")
    global openbabel_prepared
    if openbabel_prepared:
        print("openbabel already prepared")
        return

    if sys.platform != 'win32': # Only Windows seem to require this
        return

    if sys.version_info < (3, 8): # No need for this in Python 3.7 or earlier
        return

    print("We need to prepare")
    # Look for OpenBabel on the path
    paths = os.environ['PATH'].split(';') # Windows only here, no need to look for other separators
    openbabel_paths = [p for p in paths if 'openbabel' in p.lower()]
    print("All potential openbabel paths: ", openbabel_paths)
    if not openbabel_paths:
        raise ValueError("Can't locate OpenBabel in the PATH, please add it")

    # We are not picky here, just add all these directories as candidates
    for path in openbabel_paths:
        print(f"Adding {path} to DLL loading path")
        os.add_dll_directory(path)
        print("OK")
    openbabel_prepared = True
    print("Done")
    formatters.csm_out_pipe = sys.stdout