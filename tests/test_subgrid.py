# -*- coding: utf-8 -*-


import unittest
import logging

import numpy as np
import scipy.interpolate
import rasterio.transform

from flowmap import subgrid

logger = logging.getLogger(__name__)


class TestSubgrid(unittest.TestCase):

    def setUp(self):
        node_coordinates = np.array([
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
            [1.5, 0.5],
            [0.5, 1.5]
        ])
        # origami pattern for a dove
        faces = np.ma.masked_equal([
            [0, 1, 2, 3],
            [1, 4, 2, -9999],
            [2, 5, 3, -9999]
        ], -9999)
        mask = np.repeat(faces.mask[:, :, np.newaxis], repeats=2, axis=2)
        face_coordinates = node_coordinates[faces.filled(0)]
        face_coordinates = np.ma.masked_array(face_coordinates, mask)
        face_centers = np.ma.mean(face_coordinates, axis=1)

        self.grid = dict(
            node_coordinates=node_coordinates,
            faces=faces,
            face_coordinates=face_coordinates,
            face_centers=face_centers
        )
        self.dem = dict(
            band=np.ma.ones((5, 5), dtype='float32'),
            affine=rasterio.transform.Affine(
                0.5, 0.0, 0,
                0.0, -0.5, 5
            )

        )
        self.dem['world2px'] = lambda xy: np.vstack((~self.dem['affine']) * (xy[:, 0], xy[:, 1])).T.astype('int')

    def tearDown(self):
        pass

    def test_interpolate(self):
        values = np.zeros((self.grid['faces'].shape[0], 3), dtype='double')
        L = subgrid.build_interpolate(self.grid, values)
        assert isinstance(L, scipy.interpolate.LinearNDInterpolator)

    def test_build_tables(self):
        tables = subgrid.build_tables(self.grid, self.dem)
        assert len(tables) > 0


if __name__ == '__main__':
    sys.exit(unittest.main())
