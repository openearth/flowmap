FROM continuumio/miniconda3
MAINTAINER Fedor Baart <fedor.baart@deltares.nl>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# update system and install wget
RUN \
    apt-get install -y apt-utils && \
    echo "deb http://httpredir.debian.org/debian jessie-backports main non-free" >> /etc/apt/sources.list && \
    echo "deb-src http://httpredir.debian.org/debian jessie-backports main non-free" >> /etc/apt/sources.list && \
    apt-get update --fix-missing && \
    apt-get install -y ffmpeg wget unzip
# switch to python 3.5 (no gdal in 3.6)
RUN conda create -y -n py35 python=3.5 libgdal gdal jpeg=8d netcdf4 matplotlib scikit-image tqdm cython pillow click pandas
# install flowmap in the new environment
RUN /opt/conda/envs/py35/bin/pip install flowmap awscli
ENV PATH /opt/conda/bin:$PATH
# not sure what this is
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/opt/conda/envs/py35/bin/matroos_flowmap" ]
