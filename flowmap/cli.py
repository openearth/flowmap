# -*- coding: utf-8 -*-
import logging

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
    type=int
)
@click.option(
    "--dst_epsg",
    type=int
)
def generate(dataset, **kwargs):
    """Convert the dataset to a flow map."""
    klass = flowmap.formats.get_format(dataset, **kwargs)
    ds = klass(dataset, **kwargs)
    ds.animate()



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
