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
  # you can then run all commands in docker, for example
  docker run openearth/flowmap --help

You can run the software by typing the command (for now please prepend the flowmap command with `/opt/conda/envs/py35/bin`.

.. code:: bash

  # help
  flowmap --help

  # help per command
  flowmap generate --help

  # generate flowmap (for openearth/painting)
  flowmap generate delft3doutput.nc --src_epsg=28992 --dst_epsg=3857

  # export tables to nc format for faster subgrid calculations
  flowmap export --format tables --src_epsg=28992 delft3doutput.nc  aw_refi_def_asc.tiff --valid-range -10 10

  # export id grid (for faster lookups)
  flowmap export --format id_grid --src_epsg 28992 groesbeek_map.nc aw_ahn_d_asc.tiff

  # compute subgrid interpolation (for last timestep)
  flowmap subgrid delft3doutput.nc  aw_refi_def_asc.tiff --src_epsg=28992 --timestep -1

  # regrid the output to a tiff file can be done with gdal
  gdal_grid -zfield subgrid_waterlevel groesbeek_map_waterlevel_last.geojson groesbeek_map_waterlevel_last_idw.tiff -outsize 16069 20071 -a invdistnn:power=3.0:max_points=4:radius=8 -txe 188819.156 196867.156  -tye 426992.399 416956.899

  # extracting the relevant contour
  gdalwarp -q -cutline "D:/11201337 Water op Straat WS Rivierenland/008. Model/B. Results/LeerdamWest/case14/Leerdam_contour.shp" -tr 0.5 0.5 "D:/11201337 Water op Straat WS Rivierenland/008. Model/F. Post Subgrid/Leerdam/from Fedor/wd_v20180131.tif"


In the case of Delft3D you can convert the default nefis output to netCDF using the vs_trim2nc.m matlab script.
There is also direct nefis support in development, but that has not been properly tested.
By default you will want to project to the web mercator projection. Then you can reuse the velocities as pixels/s in an animation.

Exporting files
---------------
There are several export files that can be generated using the `flowmap export` command.

The `id_grid` is needed to export `tables`. The subgrid `tables` are needed for the subgrid command. The `hull` file is needed for interpolation and for flowmaps. File names are generated based on the grid name in the format: [grid_name]_[export_name].[suffix] and placed next to the grid file.



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
