===============================
Flowmap generator
===============================


.. image:: https://img.shields.io/pypi/v/flowmap.svg
        :target: https://pypi.python.org/pypi/flowmap

.. image:: https://img.shields.io/travis/openearth/flowmap.svg
        :target: https://travis-ci.org/openearth/flowmap

.. image:: https://readthedocs.org/projects/flowmap/badge/?version=latest
        :target: https://flowmap.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/openearth/flowmap/shield.svg
     :target: https://pyup.io/repos/github/openearth/flowmap/
     :alt: Updates


Command line utility to transform model output into a flowmap that can be used for games or gpu-based visualizations.


* Free software: GNU General Public License v3
* Documentation: https://flowmap.readthedocs.io.


Scripts
=======

Scripts that generate flowmaps for specific models:

- matroos_flowmap.sh (download model results for DCSM model and transform them to a flowmap)

Usage
=====

To use the software you can download the latest version using docker.
If you have docker installed you can download the software using the command:

.. code:: bash

  docker pull openearth/flowmap

You can run the software by typing the command (for now please prepend the flowmap command with `/opt/conda/envs/py35/bin`.

.. code:: bash

  # help
  docker run openearth/flowmap flowmap --help

  # help per command
  docker run openearth/flowmap flowmap generate --help

  # generate flowmap
  docker run openearth/flowmap flowmap generate delft3doutput.nc --src_epsg=28992 --dst_epsg=3857

  # compute subgrid interpolation
  docker run openearth/flowmap flowmap subgrid delft3doutput.nc  aw_refi_def_asc.tiff --src_epsg=28992

  # export tables to pickle format
  docker run openearth/flowmap flowmap export --format tables --src_epsg=28992 delft3doutput.nc  aw_refi_def_asc.tiff

In the case of Delft3D you can convert the default nefis output to netCDF using the vs_trim2nc.m matlab script.
There is also direct nefis support in development, but that has not been properly tested.
By default you will want to project to the web mercator projection. Then you can reuse the velocities as pixels/s in an animation.



Features
--------

* flowmap: animated vectorfield used for interactive particles
* streamlines: generate geojson of streamlines
* subgrid: generate subgrid waterdepth

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
