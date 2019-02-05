"""
UGrid file format support for flowmap software.
The file format is expected to follow the UGRID convention 1.0

http://ugrid-conventions.github.io/ugrid-conventions/

"""

import logging
import pickle
import pathlib

import netCDF4
import numpy as np
import pyugrid
import geojson

# used for transforming into a vtk grid and for particles
import tqdm
from tvtk.api import tvtk
import rasterio
import rasterio.mask
import rasterio.crs
import shapely.geometry

from .netcdf import NetCDF
from .. import particles
from .. import subgrid
from .. import topology
from ..dem import read_dem


logger = logging.getLogger(__name__)


class CustomEncoder(geojson.GeoJSONEncoder):
    """also convert numpy """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(CustomEncoder, self).default(obj)




class UGrid(NetCDF):
    """UGrid NetCDF format"""

    def validate(self):
        """validate a file"""
        valid = True

        with netCDF4.Dataset(self.path) as ds:
            variables = ds.variables.keys()
        # for now we hardcode the filenames. This can be replaced by the pyugrid, once released
        for var in ("mesh2d_face_nodes", "mesh2d_node_x", "mesh2d_node_x"):
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

        faces = np.ma.asanyarray(ugrid.faces)

        # Don't use these, these are circumcenters
        face_centers = ugrid.face_coordinates

        nodes = ugrid.nodes
        # should be a ragged array
        face_coordinates = np.ma.asanyarray(nodes[faces])
        face_coordinates[faces.mask] = np.ma.masked

        # recompute face centers
        face_centroids = face_coordinates.mean(axis=1)

        x = nodes[:, 0]
        y = nodes[:, 1]
        z = np.zeros_like(x)
        points = np.c_[x, y, z]
        return dict(
            face_coordinates=face_coordinates,
            face_centers=face_centers,
            face_centroids=face_centroids,
            faces=faces,
            points=points
        )

    def to_features(self):
        """convert grid to geojson features"""
        ugrid = self.ugrid
        faces = ugrid['faces']

        crs = geojson.crs.Named(
            properties={
                "name": "urn:ogc:def:crs:EPSG::{:d}".format(self.src_epsg)
            }
        )
        counts = (~faces.mask).sum(axis=1)
        face_coordinates = ugrid['face_coordinates']
        features= []
        for i, (face, count) in tqdm.tqdm(enumerate(zip(face_coordinates, counts)), desc='grid->features'):
            poly = shapely.geometry.Polygon(face[:count].tolist())
            geometry = poly.__geo_interface__
            geometry['crs'] = dict(crs)
            feature = geojson.Feature(id=i, geometry=geometry)
            features.append(feature)
        return features

    def to_polys(self):
        """convert grid to geojson features"""
        ugrid = self.ugrid
        faces = ugrid['faces']

        crs = geojson.crs.Named(
            properties={
                "name": "urn:ogc:def:crs:EPSG::{:d}".format(self.src_epsg)
            }
        )
        counts = (~faces.mask).sum(axis=1)
        face_coordinates = ugrid['face_coordinates']
        polys = []
        for i, (face, count) in tqdm.tqdm(enumerate(zip(face_coordinates, counts)), desc='grid->polys'):
            poly = shapely.geometry.Polygon(face[:count].tolist())
            polys.append(poly)
        return polys

    def to_polydata(self):
        """convert grid to a vtk polydata object"""
        ugrid = self.ugrid

        faces = ugrid['faces']
        points = ugrid['points']

        n_cells = faces.shape[0]

        cell_array = tvtk.CellArray()

        counts = (~faces.mask).sum(axis=1)
        assert faces.min() >= 0, 'expected 0 based faces'
        cell_idx = np.c_[counts, faces.filled(-999)].ravel()
        cell_idx = cell_idx[cell_idx != -999]
        cell_array.set_cells(n_cells, cell_idx)

        # fill in the properties
        polydata = tvtk.PolyData()
        polydata.points = points
        polydata.polys = cell_array
        return polydata

    def update_polydata(self, polydata, t):
        """set velocities to the polydata"""
        variables = self.velocities(t)
        ucx = variables['ucx']
        ucy = variables['ucy']
        vectors = np.c_[ucx, ucy, np.zeros_like(ucx)]
        polydata.cell_data.vectors = vectors
        polydata.cell_data.vectors.name = 'vector'
        polydata.modified()


    def waterlevel(self, t):
        """lookup the waterlevel, depth and volume on timestep t"""
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
        """lookup the velocities on the cell centers on timestep t"""
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
        """compute streamlines for timestep t"""
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
        new_name = path.with_name(path.stem + '_streamlines').with_suffix('.json')
        # save the particles
        particles.export_lines(lines, str(new_name))

    def build_is_grid(self, dem):
        """create a map in the same shape as dem denoting if a pixel is in or outside the grid"""
        is_grid = np.zeros_like(dem['band'], dtype='bool')
        polydata = self.to_polydata()
        hull = topology.concave_hull(polydata)
        # convert hull to geojson
        crs = geojson.crs.Named(
            properties={
                "name": "urn:ogc:def:crs:EPSG::{:d}".format(self.src_epsg)
            }
        )
        geometry = geojson.Polygon(coordinates=[hull.tolist()], crs=crs)
        # rasterize hull
        is_grid = rasterio.mask.geometry_mask(
            [geometry],
            out_shape=dem['band'].shape,
            transform=dem['affine'],
            # inside is grid
            invert=True
        )
        return is_grid

    def subgrid(self, t):
        """compute refined waterlevel using detailled dem, using subgrid or interpolate method"""
        dem = read_dem(self.options['dem'])
        ugrid = self.ugrid
        data = self.waterlevel(t)

        # this is slow
        table_name = self.generate_name(
            self.path,
            suffix='.nc',
            topic='tables'
        )
        # We need these two tables to do our computations
        table_path = pathlib.Path(table_name)

        # check if they exist
        if table_path.exists():
            logger.info('reading subgrid tables from %s', table_path)
            # tables are now in array format with metadata
            tables = subgrid.import_tables(str(table_path))
            metadata = tables['metadata']
            # get the valid range
            valid_range = metadata.get('valid_range')
        else:
            command = 'flowmap export --format tables {} {}'.format(
                self.path,
                self.options['dem']
            )
            msg = 'Create subgrid tables using the command: \n{}'.format(
                command
            )
            logger.warn(msg)
            raise IOError('Subgrid tables not found')

        logger.info('computing subgrid features')
        feature_collection = subgrid.compute_waterlevels(
            ugrid,
            dem,
            tables,
            data
        )
        waterlevel_name = self.generate_name(
            self.path,
            suffix='.json',
            topic='waterlevel',
            counter=t
        )
        logger.info('writing waterlevels features to {}'.format(waterlevel_name))
        # save featuress
        crs = geojson.crs.Named(
            properties={
                "name": "urn:ogc:def:crs:EPSG::{:d}".format(self.src_epsg)
            }
        )
        feature_collection['crs'] = crs
        with open(waterlevel_name, 'w') as f:
            geojson.dump(feature_collection, f, cls=CustomEncoder, allow_nan=False, ignore_nan=True)

        # interpolate waterlevels to grid (this is file based)
        interpolated_waterlevel_name = self.generate_name(
            self.path,
            suffix='.tiff',
            topic='waterlevel_idw',
            counter=t
        )
        logger.info(
            'writing interpolated waterlevels to {}'.format(
                interpolated_waterlevel_name
            )
        )
        subgrid.interpolate_waterlevels(
            waterlevel_name,
            interpolated_waterlevel_name,
            dem,
            epsg=self.src_epsg
        )
        # interpolated_waterlevel_name is now created
        with rasterio.open(str(interpolated_waterlevel_name)) as ds:
            interpolated_waterlevel = ds.read(1, masked=True)

        # buildings are now partially filtered out, values in waterlevels are interpolated into the buildings (radius 8)
        # apply invalid dem mask (filter out buildings) to leave out these cells in the waterdepth
        invalid_mask = False
        if valid_range is not None:
            invalid_mask = np.logical_or(
                dem['band'] < valid_range[0],
                dem['band'] > valid_range[1]
            )
        mask = np.logical_or(
            interpolated_waterlevel.mask,
            invalid_mask
        )

        interpolated_waterlevel = np.ma.masked_array(
            interpolated_waterlevel,
            mask
        )
        waterdepth_name = self.generate_name(
            self.path,
            suffix='.tiff',
            topic='waterdepth',
            counter=t
        )
        logger.info(
            'writing waterdepth to {}'.format(
                waterdepth_name
            )
        )
        # now for the final computation
        waterdepth = interpolated_waterlevel - dem['band']
        # save results
        subgrid.export_grid(
            waterdepth_name,
            waterdepth,
            affine=dem['affine'],
            width=dem['width'],
            height=dem['height'],
            epsg=self.src_epsg
        )

    def export(self, format, **kwargs):
        """export dataset"""
        crs = geojson.crs.Named(
            properties={
                "name": "urn:ogc:def:crs:EPSG::{:d}".format(
                    self.src_epsg
                )
            }
        )
        # TODO: export vector plot data
        if format == 'hull':
            poly = self.to_polydata()
            cells = [
                list(poly.get_cell(idx).points)
                for idx
                in range(poly.number_of_cells)
            ]
            multi_polygon = geojson.MultiPolygon(coordinates=cells, crs=crs)
            feature = geojson.Feature(id='grid', geometry=multi_polygon)
            new_name = self.generate_name(
                self.path,
                suffix='.json',
                topic=format
            )

            with open(new_name, 'w') as f:
                geojson.dump(feature, f, cls=CustomEncoder, allow_nan=False, ignore_nan=True)
        elif format == 'tables':
            dem = read_dem(self.options['dem'])

            # lookup id grid and give instructions how
            # to create it if it doesn't exist
            id_grid_name = self.generate_name(
                self.path,
                suffix='.tiff',
                topic='id_grid'
            )
            id_grid_path = pathlib.Path(id_grid_name)
            if id_grid_path.exists():
                id_grid = subgrid.import_id_grid(str(id_grid_path))
            else:
                # warning message with suggestion for import command
                command = 'flowmap export --format id_grid {} {} --src_epsg {}'.format(

                    pathlib.Path(self.path).relative_to(pathlib.Path('.').absolute()),
                    pathlib.Path(self.options['dem']).relative_to(pathlib.Path('.').absolute()),
                    self.src_epsg
                )
                msg = 'Create id_grid using the command: \n{}'.format(
                    command
                )
                logger.warn(msg)
                raise IOError('Id Grid not found')
            id_grid = subgrid.import_id_grid(id_grid_path)
            ugrid = self.ugrid
            tables = subgrid.build_tables(ugrid, dem, id_grid, kwargs.get('valid_range'))
            new_name = self.generate_name(
                self.path,
                suffix='.nc',
                topic=format
            )
            # save options to the netcdf file
            subgrid.create_export(new_name, n_cells=len(tables), n_bins=20, attributes=kwargs)
            subgrid.export_tables(new_name, tables)
        elif format == 'id_grid':
            dem = read_dem(self.options['dem'])
            polys = self.to_polys()
            id_grid = subgrid.build_id_grid(polys, dem)
            new_name = self.generate_name(
                self.path,
                suffix='.tiff',
                topic=format
            )
            subgrid.export_grid(
                new_name,
                id_grid,
                affine=dem['affine'],
                width=dem['width'],
                height=dem['height'],
                epsg=self.src_epsg
            )
        else:
            raise ValueError('unknown format: %s' % (format, ))


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
