"""
UGrid file format support for flowmap software.
The file format is expected to follow the UGRID convention 1.0

http://ugrid-conventions.github.io/ugrid-conventions/

"""

import logging
import pathlib

# TODO, switch to pyugrid after next release
import netCDF4
import numpy as np

# used for transforming into a vtk grid and for particles
from tvtk.api import tvtk

from .netcdf import NetCDF
from .. import particles

logger = logging.getLogger(__name__)


class UGrid(NetCDF):
    """UGrid NetCDF format"""

    def validate(self):
        """validate a file"""
        valid = True

        with netCDF4.Dataset(self.path) as ds:
            variables = ds.variables.keys()
        # for now we hardcode the filenames. This can be replaced by the pyugrid, once released
        for var in ("mesh2d_ucx", "mesh2d_ucy", "mesh2d_face_nodes", "mesh2d_node_x", "mesh2d_node_x"):
            if var not in variables:
                logger.warn(
                    "%s not found in variables of file %s",
                    var,
                    self.path
                )
                valid = False
                return valid
        return valid

    @property
    def grid(self):
        """return the grid variables including coordinates in space and time"""
        # hard coded names for now
        with netCDF4.Dataset(self.path) as ds:
            x = ds.variables['mesh2d_node_x'][:]
            y = ds.variables['mesh2d_node_y'][:]
            faces = ds.variables['mesh2d_face_nodes'][:]
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
        # collect all variables
        grid = dict(
            x=x,
            y=y,
            z=z,
            points=points,
            faces=faces
        )
        return grid

    def to_polydata(self):
        grid = self.grid

        faces = grid['faces']
        points = grid['points']

        n_cells = faces.shape[0]
        points.shape, n_cells

        cell_array = tvtk.CellArray()

        # for now we assume a grid with only quads.
        # switch to pyugrid to support flexible meshes (with quads + triangles)
        assert not hasattr(faces, 'mask'), 'should not be a masked array'
        assert faces.shape[1] == 4, 'we expect quads only'

        # For unstructured grids you need to count the number of edges per cell
        counts = np.ones((n_cells), dtype=faces.dtype) * 4
        cell_idx = np.c_[counts, faces - 1].ravel()
        cell_array.set_cells(n_cells, cell_idx)

        # fill in the properties
        polydata = tvtk.PolyData()
        polydata.points = points
        polydata.polys = cell_array
        return polydata

    def update_polydata(self, polydata, t):
        variables = self.variables(t)
        ucx = variables['ucx']
        ucy = variables['ucy']
        vectors = np.c_[ucx, ucy, np.zeros_like(ucx)]
        polydata.cell_data.vectors = vectors
        polydata.cell_data.vectors.name = 'vector'
        polydata.modified()

    def variables(self, t):
        with netCDF4.Dataset(self.path) as ds:
            # cumulative velocities
            ucx = ds.variables['mesh2d_ucx'][t]
            ucy = ds.variables['mesh2d_ucy'][t]
        return dict(
            ucx=ucx,
            ucy=ucy
        )

    def streamlines(self, t):
        polydata = self.to_polydata()
        self.update_polydata(polydata, t)
        seed = particles.make_particles(polydata, n=self.options.get('n_particles', 1000))
        tracer = particles.make_tracer_pipeline(polydata, seed)
        tracer.update()
        lines = particles.extract_lines(tracer)

        # convert coordinates
        def line2lonlat(line):
            lonlatz = self.srs['src2wgs84'].TransformPoints(line)
            return np.array(lonlatz)
        # replace lines with lonlat
        lines['line'] = lines['line'].apply(
            line2lonlat
        )
        # create a new name
        path = pathlib.Path(self.path)
        new_name = path.with_name(path.stem + '_streamlines').with_suffix('.geojson')
        # save the particles
        particles.export_lines(lines, str(new_name))
