# -*- coding: utf-8 -*-


import unittest

import numpy as np
from tvtk.api import tvtk

from flowmap import particles


class TestParticles(unittest.TestCase):

    def setUp(self):
        x = np.arange(1, 10)
        y = np.sqrt(x)
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
        polydata = tvtk.PolyData()
        polydata.points = points
        self.polydata = polydata


    def tearDown(self):
        pass

    def test_make_seed(self):
        seed = particles.make_particles(self.polydata)
        assert isinstance(seed, tvtk.PolyData)

    def test_make_tracer(self):
        tracer = particles.make_tracer_pipeline(
            polydata=self.polydata,
            seed=self.polydata
        )
        assert isinstance(tracer, tvtk.StreamTracer)

if __name__ == '__main__':
    sys.exit(unittest.main())
