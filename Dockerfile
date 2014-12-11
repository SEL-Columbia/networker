FROM ubuntu:latest

MAINTAINER Brandon Ogle

RUN \
  apt-get update && apt-get install -y -q \ 
  build-essential \
  make \
  wget \
  gcc \
  git 

RUN \
  wget http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh -O miniconda.sh && \
  /bin/bash miniconda.sh -b -p $HOME/anaconda/
 
ENV PATH /root/anaconda/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN \
  conda install --yes --force numpy scipy pandas networkx cython nose decorator numba dateutil shapely pytz && \
  conda install --yes -c https://conda.binstar.org/osgeo gdal && \
  conda install --yes -c https://conda.binstar.org/dougal libspatialindex && \  
  conda install --yes -c https://conda.binstar.org/dougal rtree && \
  conda install --yes -c https://conda.binstar.org/anaconda llvm && \
  conda install --yes -c https://conda.binstar.org/anaconda llvmpy && \
  mkdir code && cd code && git clone https://github.com/SEL-Columbia/NetworkBuild.git && \
  git clone https://github.com/SEL-Columbia/networkplanner.git && \
  mv NetworkBuild/networkbuild ~/anaconda/lib/python2.7/site-packages && \ 
  mv networkplanner/np ~/anaconda/lib/python2.7/site-packages && \
  cd ~/anaconda/lib/python2.7/site-packages/networkbuild && nosetests
