python setup.py prepare
python setup.py build_ext
PUSHD
ECHO copying pyd files
cd build
   for /r %%a in (*.pyd) do (
     COPY "%%a" "%~dp0%%~nxa"
   )
POPD