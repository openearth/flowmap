import logging
import pathlib

# TODO, switch to pyugrid after next release
import netCDF4
import numpy as np
from tvtk.api import tvtk
from tvtk.common import configure_input, configure_source_data, configure_input_data
import vtk

from .netcdf import NetCDF


logger = logging.getLogger(__name__)


class UGrid(NetCDF):
    """UGrid NetCDF format"""

    def validate(self):
        """validate a file"""
        valid = True

        with netCDF4.Dataset(self.path) as ds:
            variables = ds.variables.keys()
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
        with netCDF4.Dataset(self.path) as ds:
            x = ds.variables['mesh2d_node_x'][:]
            y = ds.variables['mesh2d_node_y'][:]
            faces = ds.variables['mesh2d_face_nodes'][:]
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
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

    @staticmethod
    def make_particles(polydata, n=100, mode='random'):
        points = polydata.points.to_array()
        xmax, ymax, zmax = points.max(axis=0)
        xmin, ymin, zmin = points.min(axis=0)
        # take the srt of n
        n = np.round(np.sqrt(n)).astype('int')
        if mode != 'random':
            seed_X, seed_Y = np.meshgrid(
                np.linspace(xmin, xmax, num=n),
                np.linspace(ymin, ymax, num=n)
            )
        else:
            seed_X = np.random.random((n, n)) * (xmax - xmin) + xmin
            seed_Y = np.random.random((n, n)) * (ymax - ymin) + ymin
        seed_Z = np.zeros_like(seed_X)
        seed_points = np.c_[seed_X.ravel(), seed_Y.ravel(), seed_Z.ravel()]
        seed = tvtk.PolyData()
        seed.points = seed_points
        return seed

    @staticmethod
    def make_tracer_pipeline(polydata, seed, dt=None, l=1000):
        # create elements of the pipeline
        tracer = tvtk.StreamTracer()


        # # # You can compute up to 1km per timestep (8km/hr for )

        n = 200  # max number of steps
        if dt is not None:
            # maximum velocity
            uvw = polydata.cell_data.get_array(0).to_array()
            max_velocity = np.power(uvw, 2).sum(axis=1).max()
            # s * m / s
            l = dt * max_velocity


        # maximum 1km
        tracer.maximum_propagation = l
        # # # In m
        tracer.integration_step_unit = vtk.vtkStreamTracer.LENGTH_UNIT
        # # # Minimum 5 per step
        tracer.minimum_integration_step = (l/n)
        # # # Maximum 100m per step
        tracer.maximum_integration_step = 10*(l/n)
        # # # Maximum 200 steps
        tracer.maximum_number_of_steps = n
        # # # Maximum error 1cm
        tracer.maximum_error = 1e-2
        # # # We use a path integration. You could argue that you need a
        # # # particle tracking algorithm that matches the numerical grid
        # # # (in our case edge velocities
        # # # and integration over a cell instead of over a line)
        tracer.integrator_type = 'runge_kutta45'

        tracer.debug = True

        cell2point = tvtk.CellDataToPointData()
        # setup the pipeline
        configure_input(cell2point, polydata)
        configure_input(tracer, cell2point)
        configure_source_data(tracer, seed)
        return tracer

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
        seed = self.make_particles(polydata)
        tracer = self.make_tracer_pipeline(polydata, seed)

    def meta(self):
        metadata = super().meta()
        metadata['metadata']['format'] = str(self)
        return metadata
