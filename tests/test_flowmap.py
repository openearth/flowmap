#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_flowmap
----------------------------------

Tests for `flowmap` module.
"""


import sys
import unittest
from click.testing import CliRunner

from flowmap import flowmap
from flowmap import cli


class TestFlowmap(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_flowmap(self):
        result = flowmap()
        assert result is not None

    def test_command_line_interface(self):
        runner = CliRunner()
        result = runner.invoke(cli.generate)
        assert result.exit_code == 0
        assert 'flowmap.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output


if __name__ == '__main__':
    sys.exit(unittest.main())
