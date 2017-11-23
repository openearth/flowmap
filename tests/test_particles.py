# -*- coding: utf-8 -*-


import unittest

import numpy as np
from tvtk.api import tvtk

from flowmap import particles


class TestParticles(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_make_seed(self):
        x = np.arange(1, 10)
        y = np.sqrt(x)
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
        polydata = tvtk.PolyData()
        polydata.points = points
        seed = particles.make_particles(polydata)
        assert isinstance(seed, tvtk.PolyData)


if __name__ == '__main__':
    sys.exit(unittest.main())
