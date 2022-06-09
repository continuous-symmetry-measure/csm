import platform
import sys
import os

# Install the requirements


def get_openbabel_dependency():
    plat = platform.system()
    if plat == 'Windows':
        if sys.version_info.minor == 7:
            return 'openbael==3.1.1'
        elif sys.version_info.minor == 8:
            return "openbabel-wheels/openbabel-3.1.1-cp38-cp38-win_amd64.whl"
        elif sys.version_info.minor == 9:
            return "openbabel-wheels/openbabel-3.1.1-cp39-cp39-win_amd64.whl"
        else: #elif sys.version_info.minor == 10:
            return "openbabel-wheels/openbabel-3.1.1-cp310-cp310-win_amd64.whl"

    else: # elif plat == 'Linux':
        return 'openbabel==3.1.1.1'


def get_openbabel_url():
    # Returns the the openbabel installation URL
    pass

def main():
    # Install the proper openbabel on Windows, then install requirements.windows.txt on windows.
    # On linux, just install requirements.linux.txt,
    # use os.system('pip install...') to install something.

    plat = platform.system()
    if plat == 'Windows':
        os.system('pip install ' + get_openbabel_dependency())
        os.system('pip install -r src/requirements.windows.txt' )
    
    else: # elif plat == 'Linux':
        os.system('pip install -r src/requirements.linux.txt' )

    pass

if __name__ == '__main__':
    main()