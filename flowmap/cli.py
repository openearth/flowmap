# -*- coding: utf-8 -*-
import logging

import click

logger = logging.getLogger(__name__)

@click.group()
def cli():
    pass

@cli.command()
@click.argument('dataset', type=click.Path(exists=True))
def generate(dataset):
    """Convert the dataset to a flow map."""
    logger.info("Converting dataset to flowmap. ")

@cli.command()
def formats():
    """List the available formats"""



if __name__ == "__main__":
    generate()
