# CSM Version Creation

CSM realeses new versions from time to time.  
This document explains the steps taken for a new version to be released. 

---

## Tests

- ### Feature Testing

    Each new feature requires a unit test to check the feature.  
    The tests are implemented as **pytests**.  

    First, the feature will be tested by the devloper on VSCode, and later on as a regression test on CI. 

- ### Regression Tests

    Each new featur test mentioned above will become part of the CI regression tests.
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

The build is automated by GitHub actions and triggered on a push of a new tag
This happen due to the `wheels.yml` workflow. 
The workflow is composed of a few jobs as following:
- ### Creation of the Linux wheels
    The Job runs on an Ubuntu virtual machine and uses a Docker container to build the wheels.  
    - First, the Docker image is build using the Dockerfile.manylinux file.
      The image is based on the `quay.io/pypa/manylinux_2_24_x86_64` image, and produces the wheels by executing the follwoing steps:
      - #### Installation of basic tools  
        See the Dockerfil.naylinux for details
      - #### OpenBabel installation:
        The required installation file is stored in GitHub. The image downloads and installs it.
      - #### Installation of the requirments for all python versions
      - #### Creation of the wheels
    - Next, a container is created and run to create the wheels as mentioned above.
    - And finally, the whhels are copie from the container an duploded to Release.
        
- ### Creation of the Windows wheels
    The Job runs on a Windows virtual machine and executes the following steps:
    - Checkout the lates version
    - Setup the python environment and install requirments.
      This step uses the install_requirments.py script, which determines (according to the platform and python version) which OpenBabel installation should be installed.  
      All OpenBabel installation files are stored in GitHub and installed on the virtual machine. 
    - Creation of a wheel for each python version.
    - Upload wheels to Release

- ### Creation of the source package
    This Job runs similarly to the previous Job, with the difference of the tar.gz source package created and uploaded to Release rather than a wheel for each python version.

- ### Building the Docker image
    After the wheels are built, a Docker image is created and pushed to DockerHub.
    This Job runs on an Ubuntu virtual machine and executes the following steps:
    - Checkout the lates version to get the Dockerfile
    - Set environment variables to be used later on.
    - Fetch the latest linux wheels
    - Login to DockerHub using the username and password from github secrets
    - Creation of the docker image using the Dockerfile.
      The Docker image is based on the teamcsm/openbabel:py3.9-openbabel3.1.1 image and includes the installation of the CSM wheels created on previous steps.
    - Uploading the ne Docker image to Docker Hub

> All steps are executed automatically once a new tag is added.  After completion, all the wheels can be downloaded from GitHub, and the Docker image is uploaded to DockerHub and ready to use.  

