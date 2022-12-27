# How to run locally on Windows
- install the openbabel gui http://openbabel.org/wiki/Category:Installation

- clone the CSM repo

- go to python directory, create virtual env and install requirements.txt

- run src\rebuild.bat

- run python and check that openbabel is installed correctly.  
- `python .\src\csm\main\csm_run.py --help`  

# How to run latest version with docker
`docker run -i -t -d -v testfiles --name csm_latest teamcsm/csm bin/bash`
`docker cp C:\temp\test-approx.sdf csm_latest:testfiles`
`docker attach csm_latest`