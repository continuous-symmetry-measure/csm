python setup.py prepare
python setup.py build_ext
ECHO copying pyd files
for /R .\build %%f in (*.pyd) do copy %%f .\csm