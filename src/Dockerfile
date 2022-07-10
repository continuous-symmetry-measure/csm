# A Dockerfile for CSM, with Python 3.5, Openbabel and csm installed.
# Adapted from various Dockerfiles for Openbabel and Python

# TODO: Move the openbabel 3.1.1 on Python 3.7
FROM teamcsm/openbabel:py3.7-openbabel3.1.1
LABEL by="The Research Software Company"
LABEL website="https://www.chelem.co.il"

RUN pip install numpy
RUN pip install scipy
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.5/eigen-3.3.5.tar.gz && \
	tar -xvf eigen-3.3.5.tar.gz && \
	mv /eigen-3.3.5 /usr/include/eigen335
ENV CPATH=/usr/include/eigen335

ARG WHL
ARG CSM_VER

# Switch to csm 1.2.6
RUN pip install $WHL csm==$CSM_VER
RUN pip install --upgrade numpy

