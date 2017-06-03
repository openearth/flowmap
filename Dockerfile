FROM continuumio/anaconda3
MAINTAINER Fedor Baart <fedor.baart@deltares.nl>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# update system and install wget
RUN \
    apt-get install -y apt-utils && \
    echo "deb http://httpredir.debian.org/debian jessie-backports main non-free" >> /etc/apt/sources.list && \
    echo "deb-src http://httpredir.debian.org/debian jessie-backports main non-free" >> /etc/apt/sources.list && \
    apt-get update --fix-missing && \
    apt-get install -y ffmpeg wget
# switch to python 3.5 (no gdal in 3.6)
RUN conda install -y python=3.5
RUN conda install -y gdal jpeg=8d
RUN pip install flowmap
ENV PATH /opt/conda/bin:$PATH
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
