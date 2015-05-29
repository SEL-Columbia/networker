FROM ubuntu:14.04
MAINTAINER Chris Natali

RUN \
  apt-get update && apt-get install -y -q \ 
  build-essential \
  make \
  wget \
  gcc \
  git 

RUN \
  wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh && \
  /bin/bash miniconda.sh -b -p $HOME/anaconda/
 
ENV PATH /root/anaconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN \
  conda config --add channels 'ioos' && \
  conda install --yes -c sel networkplanner-metrics && \
  conda install --yes -c sel networker && \
  ldconfig /root/anaconda/lib && \
  mkdir src && cd src && git clone https://github.com/SEL-Columbia/networker.git && \
  cd networker && nosetests
