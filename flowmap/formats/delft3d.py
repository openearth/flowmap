import logging
import functools
import itertools
import json

import pandas
import geojson
import tqdm
import netCDF4
import scipy.interpolate
import numpy as np
import skimage.draw
import skimage.morphology
import matplotlib.pyplot as plt
import matplotlib.colors

# we need a separate transform that keeps masks
from .formats import transform
from .netcdf import NetCDF
from .json_encoder import RegistryEncoder as Encoder

matplotlib.use('Agg')

logger = logging.getLogger(__name__)


class Delft3DMatlab(NetCDF):
    """NetCDF converted with vs_nefis2nc"""
    @property
    @functools.lru_cache()
    def grid(self):
        """generate grid variables"""
        logger.debug("Computing grid variables")
        src2wgs84 = self.srs['src2wgs84']
        src2utm = self.srs['src2utm']
        src2web = self.srs['src2web']
        utm2web = self.srs['utm2web']

        grid = {}
        with netCDF4.Dataset(self.path) as ds:
            time = netCDF4.num2date(ds.variables['time'][:],
                                    ds.variables['time'].units)

            x, y = ds.variables['x'][:], ds.variables['y'][:]
            # same thing for contours
            xcc = ds.variables['grid_x'][:]
            ycc = ds.variables['grid_y'][:]

            # initial values (used to determine shapes and stuff, maybe remove if not used)
            grid["s1_0"] = np.squeeze(ds.variables['waterlevel'][0])
            grid["u1_0"] = np.squeeze(ds.variables['velocity_x'][0])
            grid["v1_0"] = np.squeeze(ds.variables['velocity_y'][0])

        # store variables in grid
        grid['time'] = time
        grid['x'] = x
        grid['y'] = y

        # compute other coordinate systems
        Lon, Lat = transform(x, y, src2wgs84)
        X_web, Y_web = transform(x, y, src2web)
        X_utm, Y_utm = transform(x, y, src2utm)

        # return
        grid["lon"] = Lon
        grid["lat"] = Lat
        grid["X_web"] = X_web
        grid["Y_web"] = Y_web
        grid["X_utm"] = X_utm
        grid["Y_utm"] = Y_utm

        grid["x_src_c"] = xcc
        grid["y_src_c"] = ycc

        # is any of the contour points masked?
        mask = np.logical_or(xcc.mask, ycc.mask).any(axis=-1)
        xycc = np.stack([xcc[~mask], ycc[~mask]], axis=-1)
        # shape of the non masked contours
        old_shape = xcc[~mask].shape
        # shape of the points
        new_shape = np.prod(xycc.shape[:2]), +  xycc.shape[2],
        xy = xycc.reshape(new_shape)
        lon_c, lat_c, _ = np.asarray(src2wgs84.TransformPoints(xy)).T
        x_utm_c, y_utm_c, _ = np.asarray(src2utm.TransformPoints(xy)).T
        # reshape back
        lon_c = lon_c.reshape(old_shape)
        lat_c = lat_c.reshape(old_shape)
        x_utm_c = x_utm_c.reshape(old_shape)
        y_utm_c = y_utm_c.reshape(old_shape)
        # store in the proper locations
        Lon_c = np.ma.masked_all_like(xcc)
        Lat_c = np.ma.masked_all_like(ycc)
        Lon_c[~mask] = lon_c
        Lat_c[~mask] = lat_c
        X_utm_c = np.ma.masked_all_like(xcc)
        Y_utm_c = np.ma.masked_all_like(ycc)
        X_utm_c[~mask] = x_utm_c
        Y_utm_c[~mask] = y_utm_c
        # return
        grid["lon_c"] = lon_c
        grid["Lon_c"] = Lon_c
        grid["Lat_c"] = Lat_c
        grid["X_utm_c"] = X_utm_c
        grid["Y_utm_c"] = Y_utm_c

        grid['X_distort'], grid['Y_distort'] = self.compute_distortion(X_utm, Y_utm, utm2web)
        logger.debug("Grid generated")
        return grid

    def animate(self):
        """generate an animation (set of png's)"""

        logger.debug("Generating animation")
        F, X, Y = self.canvas['F'], self.canvas['X'], self.canvas['Y']

        X_web = self.grid['X_web']
        Y_web = self.grid['Y_web']
        coordinate_mask = np.logical_or(X_web.mask, Y_web.mask)

        framescale = float(self.framescale)
        count = itertools.count()

        N = matplotlib.colors.Normalize(self.vmin, self.vmax, clip=True)
        for i in tqdm.tqdm(range(self.grid['time'].shape[0] - 1)):
            # get data for t = t
            data_0 = self.variables(i)
            data_1 = self.variables(i+1)
            u0, v0 = data_0['u1'], data_0['v1']
            u1, v1 = data_1['u1'], data_1['v1']
            for j in tqdm.tqdm(range(int(framescale))):
                u = (1.0 - (j/framescale)) * u0 + (j/framescale) * u1
                v = (1.0 - (j/framescale)) * v0 + (j/framescale) * v1
                # compute uv and mask with coordinate mask
                uv = np.c_[u[~coordinate_mask], v[~coordinate_mask].ravel()]
                F.values = uv.astype(F.values.dtype)
                UV = F(X, Y)
                RG = N(UV)
                R, G = RG[..., 0], RG[..., 1]

                # cells without a velocity
                value_mask = np.logical_and(UV[..., 0] == 0.0, UV[..., 1] == 0.0)
                # masked cells
                B = np.zeros_like(R) + np.logical_or(~self.canvas['is_grid'], value_mask)
                RGB = np.dstack([R, G, B])
                # store in filename
                # TODO: generate with ffmpeg
                plt.imsave('test_%06d.png' % (next(count),), RGB)
        logger.debug("Animation generated")

    def variables(self, t):
        with netCDF4.Dataset(self.path) as ds:
            u1 = np.squeeze(ds.variables['velocity_x'][t])
            v1 = np.squeeze(ds.variables['velocity_y'][t])
        return dict(
            u1=u1,
            v1=v1
        )

    @property
    @functools.lru_cache()
    def canvas(self):
        """determine the rendering canvas and compute coordinates"""
        logger.debug("Computing canvas properties")
        web2wgs84 = self.srs['web2wgs84']
        utm2web = self.srs['utm2web']

        grid = self.grid

        ll_web = grid['X_web'].min(), grid['Y_web'].min()
        ur_web = grid['X_web'].max(), grid['Y_web'].max()

        # we want a big map as the base layer.
        nx = 1024
        ny = 1024

        x_web_canvas = np.linspace(ll_web[0], ur_web[0], num=nx)
        y_web_canvas = np.linspace(ll_web[1], ur_web[1], num=ny)

        ll_wgs84 = web2wgs84.TransformPoint(x_web_canvas[0], y_web_canvas[0])
        ur_wgs84 = web2wgs84.TransformPoint(x_web_canvas[-1], y_web_canvas[-1])

        # coordinates in map space
        X_web_canvas, Y_web_canvas = np.meshgrid(x_web_canvas, y_web_canvas)

        X_web_c, Y_web_c = transform(grid['X_utm_c'], grid['Y_utm_c'], utm2web)

        mask = np.logical_or(X_web_c.mask.any(axis=-1), Y_web_c.mask.any(axis=-1))
        old_shape = X_web_c.shape
        new_shape = np.prod(old_shape[:2]), + old_shape[2],

        x_web_c = X_web_c.reshape(new_shape)[~mask.flatten()]
        y_web_c = Y_web_c.reshape(new_shape)[~mask.flatten()]

        x_px_c = nx * (x_web_c - ll_web[0]) / (ur_web[0] - ll_web[0])
        y_px_c = ny * (y_web_c - ll_web[1]) / (ur_web[1] - ll_web[1])

        is_grid = np.zeros((ny, nx), dtype='bool')
        for x_, y_ in zip(x_px_c, y_px_c):
            rr, cc = skimage.draw.polygon(y_, x_, is_grid.shape)
            is_grid[rr, cc] = True
        # want to dilate the grid a bit so colors will run through
        # is_grid = ~skimage.morphology.dilation(is_grid, skimage.morphology.square(5))

        # compute interpolation function
        X_web = grid['X_web']
        Y_web = grid['Y_web']
        mask = np.logical_or(X_web.mask, Y_web.mask)
        XY_web = np.c_[X_web[~mask], Y_web[~mask]]
        uv = np.c_[grid['u1_0'][~mask], grid['v1_0'][~mask]]
        F = scipy.interpolate.LinearNDInterpolator(
            XY_web,
            uv
        )
        F.fill_value = 0.0
        logger.debug("Canvas generated")
        return dict(
            X=X_web_canvas,
            Y=Y_web_canvas,
            nx=nx,
            ny=ny,
            bbox_wgs84=ll_wgs84 + ur_wgs84,
            bbox_web=ll_web + ur_web,
            is_grid=is_grid,
            mask=mask,
            F=F
        )

    def timeseries(self, i, j):
        with netCDF4.Dataset(self.path) as ds:
            u1 = np.squeeze(ds.variables['velocity_x'][..., i, j])
            v1 = np.squeeze(ds.variables['velocity_y'][..., i, j])
            s1 = np.squeeze(ds.variables['waterlevel'][:, i, j])

        date = [str(x) for x in self.grid['time']]

        df = pandas.DataFrame(
            data=dict(
                date=date,
                t=self.grid['time'],
                u1=u1,
                v1=v1,
                s1=s1
            )
        )
        return df

    def extract_points(self, points, filename="points.json"):
        grid = self.grid
        features = []
        for p in points:
            lat_i, lon_i = p
            distance = np.sqrt((grid["lat"] - lat_i)**2 + (grid["lon"] - lon_i)**2)
            i = np.argmin(distance)
            i, j = np.unravel_index(i, distance.shape)
            logger.info("distance %s", distance.min())
            logger.info("closest point for %s is %s, %s", p, i, j)

            point = geojson.Point(coordinates=[float(lon_i), float(lat_i)])

            ts = self.timeseries(i, j)
            # convert forth and back to json
            feature = geojson.Feature(id="{}_{}".format(i, j), geometry=point, properties={"series": ts})
            features.append(feature)
        feature_collection = geojson.FeatureCollection(features)
        with open(filename, "w") as f:
            json.dump(feature_collection, f, cls=Encoder)

    @staticmethod
    def compute_distortion(x_utm, y_utm, utm2web):
        """compute the advection distortion due to reprojection"""
        # compute the local distortion
        # how many meters is 1 meter in the projected system, given x,y
        x_plus_half = transform(x_utm + 0.5, y_utm, utm2web)[0]
        x_minus_half = transform(x_utm - 0.5, y_utm, utm2web)[0]
        y_plus_half = transform(x_utm, y_utm + 0.5, utm2web)[1]
        y_minus_half = transform(x_utm, y_utm - 0.5, utm2web)[1]

        # compute the deformation factor
        x_distort = (x_plus_half - x_minus_half)
        y_distort = (y_plus_half - y_minus_half)
        return x_distort, y_distort

    def validate(self):
        """check if variable is of expected format"""
        with netCDF4.Dataset(self.path) as ds:
            observed = set(ds.variables.keys())
        expected = {"time", "x", "y", "grid_x",
                    "grid_y", "waterlevel",
                    "velocity_x", "velocity_y"}
        if expected - observed:
            logger.warn("missing variables %s", expected - observed)
            return False
        return True
