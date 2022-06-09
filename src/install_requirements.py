import platform
import sys
import os

# Install the requirements

# The WINDOWS_WHEEL_DIR is at repo-root/openbabel-wheels
def get_openbabel_dependency():
    plat = platform.system()
    if plat == 'Windows':
        if sys.version_info.minor == 7:
            return 'openbabel==3.1.1'
        elif sys.version_info.minor == 8:
            return "openbabel@https://github.com/continuous-symmetry-measure/openbabel-wheels/releases/download/v1/openbabel-3.1.1-cp310-cp310-win_amd64.whl"
        elif sys.version_info.minor == 9:
            return "openbabel@https://github.com/continuous-symmetry-measure/openbabel-wheels/releases/download/v1/openbabel-3.1.1-cp39-cp39-win_amd64.whl"
        else: #elif sys.version_info.minor == 10:
            return "openbabel@https://github.com/continuous-symmetry-measure/openbabel-wheels/releases/download/v1/openbabel-3.1.1-cp38-cp38-win_amd64.whl"

    else: # elif plat == 'Linux':
        return 'openbabel==3.1.1.1'

def main():
    plat = platform.system()
    if plat == 'Windows':
        os.system('pip install ' + get_openbabel_dependency())
        os.system('pip install -r requirements.windows.txt' )
    
    else: # elif plat == 'Linux':
        os.system('pip install -r requirements.linux.txt' )

    pass

if __name__ == '__main__':
    main()