# -*- coding: utf-8 -*-
import logging

import click

from .formats import __all__ as available_formats

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
@click.argument(
    "output_dir",
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True
    )
)
def generate(dataset, output_dir):
    """Convert the dataset to a flow map."""
    logger.debug("Converting dataset %s to flowmap. Saving results in %s",
                 dataset, output_dir)




@cli.command()
def formats():
    """List the available formats"""
    logger.info("The following formats are available:")
    for format in available_formats:
        logger.info("%s", format)


if __name__ == "__main__":
    cli()
