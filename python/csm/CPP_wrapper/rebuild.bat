del /q build
python setup.py build_ext
move build\lib.win-amd64-3.5\*.pyd ..
