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
import pyugrid

# used for transforming into a vtk grid and for particles
import tqdm
import skimage.draw
from tvtk.api import tvtk
import rasterio
import rasterio.crs

from .netcdf import NetCDF
from .. import particles
from .. import subgrid
from ..dem import read_dem


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
    def ugrid(self):
        """Generate a ugrid grid from the input"""
        # TODO, lookup mesh name
        ugrid = pyugrid.UGrid.from_ncfile(self.path, 'mesh2d')
        faces = ugrid.faces
        faces_masked = np.ma.masked_array(faces, mask=not faces.fill)
        face_centers = ugrid.face_coordinates
        nodes = ugrid.nodes
        face_coordinates = nodes[faces]
        x = nodes[:, 0]
        y = nodes[:, 1]
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
        return dict(
            face_coordinates=face_coordinates,
            face_centers=face_centers,
            faces=faces_masked,
            points=points
        )

    def to_polydata(self):
        """convert grid to polydata"""
        grid = self.ugrid

        faces = grid['faces']
        points = grid['points']

        n_cells = faces.shape[0]
        points.shape, n_cells

        cell_array = tvtk.CellArray()

        counts = (~faces.mask).sum(axis=1)
        faces_zero_based = faces - 1
        cell_idx = np.c_[counts, faces_zero_based.filled(-999)].ravel()
        cell_idx = cell_idx[cell_idx != -999]
        cell_array.set_cells(n_cells, cell_idx)

        # fill in the properties
        polydata = tvtk.PolyData()
        polydata.points = points
        polydata.polys = cell_array
        return polydata

    def update_polydata(self, polydata, t):
        variables = self.velocities(t)
        ucx = variables['ucx']
        ucy = variables['ucy']
        vectors = np.c_[ucx, ucy, np.zeros_like(ucx)]
        polydata.cell_data.vectors = vectors
        polydata.cell_data.vectors.name = 'vector'
        polydata.modified()


    def waterlevel(self, t):
        # TODO: inspect mesh variable
        with netCDF4.Dataset(self.path) as ds:
            s1 = ds.variables['mesh2d_s1'][t]
            waterdepth = ds.variables['mesh2d_waterdepth'][t]
            vol1 = ds.variables['mesh2d_vol1'][t]
        return dict(
            s1=s1,
            vol1=vol1,
            waterdepth=waterdepth
        )

    def velocities(self, t):
        # TODO: inspect mesh variables
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

    def build_is_grid(self, raster):
        counts = (~self.ugrid['faces'].mask).sum(axis=1)
        is_grid = np.zeros_like(raster['band'], dtype='bool')
        polys = np.array([raster['world2px'](xy) for xy in self.ugrid['face_coordinates']])
        for i, poly in tqdm.tqdm(enumerate(polys)):
            # drawing grid mask
            # TODO: check for triangles here...
            rr, cc = skimage.draw.polygon(poly[:counts[i], 1], poly[:counts[i], 0])
            is_grid[rr, cc] = True
        return is_grid


    def subgrid(self, t, method):
        """compute refined waterlevel using detailled dem, using subgrid or interpolate method"""
        dem = read_dem(self.options['dem'])
        grid = self.ugrid
        data = self.waterlevel(t)
        if method == 'subgrid':

            logger.info('creating subgrid tables')
            # this is slow
            tables = subgrid.build_tables(grid, dem)
            logger.info('computing subgrid band')
            # this is also slow
            band = subgrid.compute_band(grid, dem, tables, data)
        elif method == 'interpolate':
            values = np.c_[data['s1'], data['vol1'], data['waterdepth']]
            # create a grid mask for the dem
            is_grid = self.build_is_grid(dem)
            logger.info('building interpolation')
            L = subgrid.build_interpolate(grid, values)
            logger.info('computing interpolation for current timestep')
            interpolated = subgrid.compute_interpolated(L, dem, data)
            band = interpolated['masked_waterdepth']
            # mask out non grid pixels
            band.mask = np.logical_or(band.mask, ~is_grid)
        else:
            raise ValueError('unknown method')
        logger.info('writing subgrid band')
        # use extreme value as nodata
        nodata = np.finfo(band.dtype).min
        options = dict(
            dtype=str(band.dtype),
            nodata=nodata,
            count=1,
            compress='lzw',
            tiled=True,
            blockxsize=256,
            blockysize=256,
            driver='GTiff',
            affine=dem['affine'],
            width=dem['width'],
            height=dem['height'],
            crs=rasterio.crs.CRS({'init': 'epsg:%d' % (self.src_epsg)})
        )
        new_name = self.generate_name(
            self.path,
            suffix='.tiff',
            topic=method,
            counter=t
        )
        with rasterio.open(str(new_name), 'w', **options) as dst:
            dst.write(band.filled(nodata), 1)

    def export(self, format):
        grid = self.ugrid



    @staticmethod
    def generate_name(path, suffix, topic=None, counter=None):
        path = pathlib.Path(path)
        base = path.stem
        if topic:
            base += '_' + topic
        if counter is not None:
            if counter == -1:
                base += '_last'
            else:
                base += '_%06d' % (counter, )
        new_name = path.with_name(base).with_suffix(suffix)
        return new_name
