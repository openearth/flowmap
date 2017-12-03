# -*- coding: utf-8 -*-


import unittest

import numpy as np
import rasterio.transform

from flowmap import subgrid


class TestSubgrid(unittest.TestCase):

    def setUp(self):
        node_coordinates = np.array([
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
            [1.5, 0.5]
        ])
        faces = np.ma.masked_equal([
            [0, 1, 2, 3],
            [1, 4, 2, -9999]
        ], -9999)
        mask = np.repeat(faces.mask[:, :, np.newaxis], repeats=2, axis=2)
        face_coordinates = node_coordinates[faces.filled(0)]
        face_coordinates = np.ma.masked_array(face_coordinates, mask)

        self.grid = dict(
            node_coordinates=node_coordinates,
            faces=faces,
            face_coordinates=face_coordinates
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

    def test_build_tables(self):
        tables = subgrid.build_tables(self.grid, self.dem)
        assert len(tables) > 0


if __name__ == '__main__':
    sys.exit(unittest.main())
