# -*- coding: utf-8 -*-
import logging
import json

import click

import flowmap.formats

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


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
    with open("model.json", "w") as f:
        json.dump(ds.canvas['bbox_wgs84'], f)
    ds.animate()


@cli.command()
@click.argument(
    "dataset",
    type=click.Path(
        exists=True,
        resolve_path=True
    )
)
@click.option(
    "-p",
    type=(float, float),
    multiple=True,
    required=True,
    help="point (lon, lat)"
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
