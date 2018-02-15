# -*- coding: utf-8 -*-


import unittest
import logging

import numpy as np
import rasterio.transform

from flowmap import subgrid
import flowmap.formats.ugrid

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
        face_centroids = np.ma.mean(face_coordinates, axis=1)

        ugrid = dict(
            node_coordinates=node_coordinates,
            faces=faces,
            face_coordinates=face_coordinates,
            face_centroids=face_centroids
        )

        # Mock grid
        class Grid():
            # some shared functions
            src_epsg = 28992
            to_polys = flowmap.formats.ugrid.UGrid.to_polys
        grid = Grid()
        grid.ugrid = ugrid
        self.grid = grid

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

    def test_build_id_grid(self):
        polys = self.grid.to_polys()
        id_grid = subgrid.build_id_grid(polys, self.dem)
        assert len(id_grid) > 0

    def test_build_tables(self):
        polys = self.grid.to_polys()
        id_grid = subgrid.build_id_grid(polys, self.dem)
        tables = subgrid.build_tables(self.grid.ugrid, self.dem, id_grid)
        assert len(tables) > 0

    def test_subgrid(self):
        polys = self.grid.to_polys()
        id_grid = subgrid.build_id_grid(polys, self.dem)
        tables = subgrid.build_tables(self.grid.ugrid, self.dem, id_grid)
        assert len(tables) > 0


if __name__ == '__main__':
    sys.exit(unittest.main())
