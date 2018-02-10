# -*- coding: utf-8 -*-
import logging
import json
import datetime
import pathlib

import numpy as np
import click
import uuid

import flowmap.formats


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def json_serializer(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, datetime.datetime):
        serial = obj.isoformat()
        return serial
    if isinstance(obj, np.ndarray):
        serial = obj.tolist()
        return serial
    if isinstance(obj, uuid.UUID):
        return str(obj)
    raise TypeError("Object %s of type %s not serializable " % (
        obj, type(obj))
    )


@click.group()
def cli():
    pass


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--src_epsg",
    type=int,
    required=True
)
def dump(dataset, **kwargs):
    """show some info about the dataset"""
    click.echo(kwargs)
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    click.echo(ds.dump())


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--src_epsg",
    type=int,
    required=True
)
@click.option(
    "--dst_epsg",
    type=int,
    required=True
)
@click.option(
    "--vmin",
    type=float
)
@click.option(
    "--vmax",
    type=float
)
@click.option(
    "--framescale",
    type=float,
    default=30.0
)
def generate(dataset, **kwargs):
    """Convert the dataset to a flow map."""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    ds.animate()


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
def meta(dataset, **kwargs):
    """generate metadata json file"""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    meta = ds.meta()
    path = pathlib.Path(dataset)
    path_name = str(path.with_suffix('.json'))
    with open(path_name, "w") as f:
        json.dump(meta, f, default=json_serializer)


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--src_epsg",
    type=int
)
@click.option(
    "-p",
    type=(float, float),
    multiple=True,
    required=True,
    help="point (lat, lon)"
)
def timeseries(dataset, p, **kwargs):
    """Extract timeseries from the dataset."""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    logger.info("extracting points %s", p)
    ds.extract_points(p)


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--timestep",
    type=int,
    default=-1
)
@click.option(
    "--src_epsg",
    type=int,
    required=True
)
@click.option(
    "--n_particles",
    type=int,
    default=1000
)
def streamlines(dataset, timestep, **kwargs):
    """Extract streamlines from the dataset. By default based on the last timestep"""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    logger.info("extracting streamlines")
    if hasattr(ds, 'streamlines'):
        ds.streamlines(timestep)
    else:
        raise ValueError('streamlines not yet supported for format', klass)


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.argument(
    "dem",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--timestep",
    type=int,
    default=-1
)
@click.option(
    "--method",
    type=click.Choice(['waterdepth', 'waterlevel']),
    default="waterlevel"
)
@click.option(
    "--src_epsg",
    type=int,
    required=True
)
def subgrid(dataset, dem, timestep, method, **kwargs):
    """Create a-posteriori subgrid map for the dataset. By default based on the last timestep. If method `waterdepth` is used it generates a geojson file with subgrid waterdepth per cell. If method `waterlevel` is used it also generates a tiff file with subgrid waterlevels (interpolated waterdepth - dem) per pixel."""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, dem=dem, **kwargs)
    logger.info("extracting subgrid")
    if hasattr(ds, 'subgrid'):
        ds.subgrid(timestep, method=method)
    else:
        raise ValueError('subgrid not yet supported for format', klass)

@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.argument(
    "dem",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "--format",
    type=click.Choice(['hull', 'tables', 'id_grid'])
)
@click.option(
    "--src_epsg",
    type=int,
    required=True
)
@click.option(
    "--valid-range",
    nargs=2,
    type=click.Tuple([float, float]),
    default=(None, None)
)
def export(dataset, dem, format, **kwargs):
    """Export a file. The `id_grid` is needed to export tables.
    The subgrid `tables` are needed for the subgrid command.
    The `hull` file is needed for interpolation and for flowmaps.
    File names are generated based on the grid name in the format:
    [grid_name]_[export_name].[suffix]
    """
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, dem=dem, **kwargs)
    logger.info("exporting grid")
    if hasattr(ds, 'export'):
        ds.export(format=format, **kwargs)
    else:
        raise ValueError('exporting not yet supported for format', klass)

@cli.command()
def formats():
    """List the available formats"""
    available_formats = []
    for name in dir(flowmap.formats):
        format = getattr(flowmap.formats, name)
        # looks like a format
        if hasattr(format, "dump"):
            available_formats.append(name)
    msg = "The following formats are available: {}".format(available_formats)
    logger.info(msg)


if __name__ == "__main__":
    cli()
