# -*- coding: utf-8 -*-
import logging
import json
import datetime

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
    raise TypeError("Object %s of type %s not serializable " % (obj, type(obj)) )


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
    type=int
)
@click.option(
    "--dst_epsg",
    type=int
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
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    with open("model.json", "w") as f:
        json.dump(ds.meta(), f, default=json_serializer)


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
def formats():
    """List the available formats"""
    logger.info("The following formats are available:")
    for format in flowmap.formats:
        # looks like a format
        if hasattr(format, "dump"):
            logger.info("%s", format)


if __name__ == "__main__":
    cli()
