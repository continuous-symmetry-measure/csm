# CSM Version Creation

CSM releases new versions from time to time.  
This document explains the steps taken for a new version to be released. 

## How to update version:
- Run all the test (see the [Tests](#tests) section bellow).
- Choose the new tag, you can [see the previous tags here](https://github.com/continuous-symmetry-measure/csm/tags).  
- Optional: Update the file src\csm\version.py with the new tag, without the 'v' prefix. (The github action will update the version.py file automatically on the release according to the tag, but you can update it here too for the development process)
- Commit and push all of your change.
- Run the code in your terminal from the csm code folder:  
```
git tag [the-new-tag, for example: `v1.3.4`]
git push --tag
```
- The git-action will start to build the new release, and run all the test on linux. You can [see the process here](https://github.com/continuous-symmetry-measure/csm/actions).
- Your new release can be found on [The Releases page](https://github.com/continuous-symmetry-measure/csm/releases)  


---

## Tests

### Features Testing

Each new feature requires a unit test to check the feature.  
The tests are implemented with *pytests*. 

First, the feature will be tested by the developer on VSCode or in Terminal, and later on as a regression test on CI. 
Running the pytests requires a few simple steps:

**Running tests in Terminal:**
- In the Terminal, create a virtual environment:  
`py -3.9 -m venv env --prompt "csm"`
- Activate the environment:  
  (Windows: `Activate-Virtenv`, Linux: `. env/bin/activate`)
- Once you are inside the virtual environment, run the requirements:  
  `cd src`
  `.\install_requirements.py`  
  `pip install pytest`  
- Build (build the c++ classes and copy the result .pyd files to the destination folder):  
  `.\rebuild.bat` (on linux run the commands from this file)
- Run the tests:
  `python -m pytest .\tests\`  
    
**Running tests on Visual Studio Code:**  
- Open the project directory in Visual Studio Code.
- In the Terminal, create a virtual environment:  
    `py -3.9 -m venv env --prompt "csm"`
- Activate the environment:  
    `Activate-Virtenv`
- Once you are inside the virtual environment, run the requirments:  
    `cd src`
    `.\install_requirements.py`
- Build:  
    `.\rebuild.bat`
- Configure VSCode to work with the virtual env
  `CTRL+SHIFT+P`  
  Select `Python: Select Interpeter`  
  Choose the Python version you are using with the virual env. ('env':venv)
- Configure Tests:
  `CTRL+SHIFT+P` 
  Select `Python: Configure Tests`   
  Select `pytest`  
  Select `. Root Directory`
- Discover and run the tests:
  on the left pane, click on the test button and refresh.
- Run the required tests/all tests.



### Regression Tests

Each new feature test mentioned above will become part of the CI regression tests.
On every new *push*, all tests will run on GitHub as a GitHub action.

This happens due to the `tests.yml` workflow, triggered on every *push*.  
The workflow works on an ubuntu virtual machine and runs the following steps:
- Checkout the latest version
- Python installation
- OpenBabel Installation
- Build CSM
- Run all pytests

## Building the Wheels and Docker Image

Once the version has been tested and is ready for release, a new tag is added and pushed to GitHub.
Tagging a version triggers creation of new wheels for the release.
Several assets are created for the new version:
- For each python version from 3.7 up to 3.10 two different wheels:
    - Windows wheels
    - Linux wheels
- A tar.gz file including the source

The build is automated by GitHub actions and triggered on a push of a new tag.  
This happens due to the `wheels.yml` workflow. 
The workflow is composed of a few jobs as following:
### Creation of the Linux wheels
The Job runs on an Ubuntu virtual machine and uses a Docker container to build the wheels.  
- First, the Docker image is build using the Dockerfile.manylinux file.
  The image is based on the `quay.io/pypa/manylinux_2_24_x86_64` image, and produces the wheels by executing the following steps:
- Installation of basic tools (See the Dockerfile.manylinux for details)
- OpenBabel installation: The required installation file is stored in GitHub. The image downloads and installs it.
- Installation of the requirments for all python versions
- Creation of the wheels.
- Next, a container is created and run to create the wheels as mentioned above.
- And finally, the wheels are copied from the container and uploaded to Release.
        
### Creation of the Windows wheels
The Job runs on a Windows virtual machine and executes the following steps:
- Checkout the latest version
- Setup the python environment and install requirments.
  This step uses the install_requirments.py script, which determines (according to the platform and python version) which OpenBabel installation should be installed.  
  All relevant OpenBabel installation files are stored in GitHub and installed on the virtual machine. 
- Creation of a wheel for each python version.
- Upload wheels to Release

### Creation of the source package
This Job runs similarly to the previous Job, with the difference of the tar.gz source package created and uploaded to Release rather than a wheel for each python version.

### Building the Docker image
After the wheels are built, a Docker image is created and pushed to DockerHub.
This Job runs on an Ubuntu virtual machine and executes the following steps:
- Checkout the latest version to get the Dockerfile.
- Set environment variables to be used later on.
- Fetch the latest linux wheels.
- Login to DockerHub using the username and password from github secrets.
- Creation of the docker image using the Dockerfile.
  The Docker image is based on the teamcsm/openbabel:py3.9-openbabel3.1.1 image and includes the installation of the CSM wheels created on previous steps.
- Uploading the new Docker image to Docker Hub

> All steps are executed automatically once a new tag is added.  After completion, all the wheels can be downloaded from GitHub, and the Docker image is uploaded to DockerHub and ready to use.  

