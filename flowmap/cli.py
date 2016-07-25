# -*- coding: utf-8 -*-
import logging

import click

logger = logging.getLogger(__name__)

@click.group()
def cli():
    pass

@cli.command()
@click.argument('dataset', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True))
def generate(dataset, output_dir):
    """Convert the dataset to a flow map."""
    logger.info("Converting dataset %s to flowmap. Saving results in %s", dataset, output_dir)

@cli.command()
def formats():
    """List the available formats"""



if __name__ == "__main__":
    generate()
