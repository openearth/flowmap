FROM continuumio/miniconda3
MAINTAINER Fedor Baart <fedor.baart@deltares.nl>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# update system and install wget
RUN \
    apt-get update --fix-missing && \
    apt-get install -y ffmpeg wget unzip libglu1-mesa-dev gcc nfs-common && \
    rm -rf /var/lib/apt/lists/*

RUN conda config --add channels conda-forge
RUN conda install libgdal gdal vtk mayavi netcdf4 numpy scipy scikit-image matplotlib pandas rasterio Shapely pillow cython

# install flowmap
COPY ./ app/
RUN pip install pip --upgrade
RUN pip install -r app/requirements_dev.txt
RUN pip install -r app/requirements.txt
RUN cd app/ && \
	python setup.py install

CMD [ "/opt/conda/bin/matroos_flowmap" ]
