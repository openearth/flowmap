import logging
import itertools
import json
import functools
import pathlib

import geojson
import tqdm
import pandas
import netCDF4
import numpy as np
import matplotlib.colors
import skimage.draw
import scipy.interpolate

from .formats import transform, points2contours, contours2vertices
from .netcdf import NetCDF

import matplotlib.pyplot as plt
matplotlib.use('Agg')

logger = logging.getLogger(__name__)


class Matroos(NetCDF):
    """FEWS Matroos format"""

    s = np.s_[:, :]

    @property
    @functools.lru_cache()
    def grid(self):
        """generate global variables"""

        src2wgs84 = self.srs['src2wgs84']
        wgs842utm = self.srs['wgs842utm']
        src2utm = self.srs['src2utm']
        src2web = self.srs['src2web']
        utm2web = self.srs['utm2web']

        with netCDF4.Dataset(self.path) as ds:
            time = netCDF4.num2date(
                ds.variables['time'][:],
                ds.variables['time'].units
            )
            analysis_time = netCDF4.num2date(
                ds.variables['analysis_time'][:],
                ds.variables['analysis_time'].units
            )
            # these lat,lon variables are actually contours of the original cells
            lat = ds.variables['lat'][self.s]
            lon = ds.variables['lon'][self.s]

            # initial values (used to determine shapes and stuff, maybe remove if not used)
            sep_0 = ds.variables['sep'][0][self.s][1:, 1:]
            u1_0 = ds.variables['velu'][0][self.s][1:, 1:]
            v1_0 = ds.variables['velv'][0][self.s][1:, 1:]

        vertices = contours2vertices(lat, lon)

        logger.info("computed vertices of shape %s", vertices.shape)
        lon_c = vertices[:, :, 0]
        lat_c = vertices[:, :, 1]
        X_utm_c, Y_utm_c = transform(lon_c, lat_c, wgs842utm)

        # coordinates are cell centers
        X_src = lon_c.mean(axis=1).reshape(u1_0.shape)
        Y_src = lat_c.mean(axis=1).reshape(u1_0.shape)
        X_web, Y_web = transform(X_src, Y_src, src2web)
        X_utm, Y_utm = transform(X_src, Y_src, src2utm)

        variables = dict(
            time=time,
            analysis_time=analysis_time,
            u1_0=u1_0,
            v1_0=v1_0,
            sep_0=sep_0,
            lat=lat,
            lon=lon,
            X_web=X_web,
            Y_web=Y_web,
            X_utm=X_utm,
            Y_utm=Y_utm,
            X_utm_c=X_utm_c,
            Y_utm_c=Y_utm_c
        )
        return variables

    @property
    @functools.lru_cache()
    def canvas(self, bbox_wgs84=None):
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

        x_web_c = X_web_c[~mask]
        y_web_c = Y_web_c[~mask]

        x_px_c = nx * (x_web_c - ll_web[0]) / (ur_web[0] - ll_web[0])
        y_px_c = ny * (y_web_c - ll_web[1]) / (ur_web[1] - ll_web[1])

        is_grid = np.zeros((ny, nx), dtype='bool')
        for x_, y_ in tqdm.tqdm(zip(x_px_c.filled(), y_px_c.filled()), desc="drawing is_grid", total=x_px_c.shape[0]):
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

    def animate(self):
        logger.debug("Generating animation")
        F, X, Y = self.canvas['F'], self.canvas['X'], self.canvas['Y']

        X_web = self.grid['X_web']
        Y_web = self.grid['Y_web']
        coordinate_mask = np.logical_or(X_web.mask, Y_web.mask)

        framescale = float(self.framescale)
        count = itertools.count()

        N = matplotlib.colors.Normalize(self.vmin, self.vmax, clip=True)
        for i in tqdm.tqdm(range(self.grid['time'].shape[0] - 1), desc="an-time"):
            data_0 = self.variables(i)
            data_1 = self.variables(i+1)
            u0, v0 = data_0['u1'], data_0['v1']
            u1, v1 = data_1['u1'], data_1['v1']

            # points for which uv
            uv0 = np.c_[u0[~coordinate_mask].filled(0), v0[~coordinate_mask].filled(0).ravel()]
            F.values = uv0.astype(F.values.dtype)
            UV0 = F(X, Y)

            uv1 = np.c_[u1[~coordinate_mask].filled(0), v1[~coordinate_mask].filled(0).ravel()]
            F.values = uv1.astype(F.values.dtype)
            UV1 = F(X, Y)

            # cells without a velocity
            value_mask = np.logical_or(
                np.logical_and(
                    UV0[..., 0] == 0.0,
                    UV0[..., 1] == 0.0
                ),
                np.logical_and(
                    UV1[..., 0] == 0.0,
                    UV1[..., 1] == 0.0
                )
            )
            mask = np.logical_or.reduce([
                ~self.canvas['is_grid']
            ])
            # value_mask

            for j in tqdm.tqdm(range(int(framescale))):
                UV = (1.0 - (j/framescale)) * UV0 + (j/framescale) * UV1
                # set U and V to 0
                UV[mask, 0] = 0.0
                UV[mask, 1] = 0.0

                # compute uv and mask with coordinate mask
                RG = N(UV)
                R, G = RG[..., 0], RG[..., 1]

                B = np.zeros_like(R) + mask
                RGB = np.dstack([R, G, B])
                # store in filename
                # TODO: generate with ffmpeg
                path = pathlib.Path(self.path)
                filename = path.parent / (path.stem + '_%06d.png' % (next(count), ))
                plt.imsave(str(filename), RGB)

        return

    def extract_points(self, points, filename="timeseries.json"):
        grid = self.grid
        records = []
        for p in points:
            lat_i, lon_i = p
            distance = np.sqrt((grid["lat"] - lat_i)**2 + (grid["lon"] - lon_i)**2)
            i = np.argmin(distance)
            i, j = np.unravel_index(i, distance.shape)
            logger.info("distance %s", distance.shape)
            logger.info("closest point for %s is %s, %s", p, i, j)
            ts = self.timeseries(i, j)

            # convert forth and back to json
            data = json.loads(ts.to_json(orient="records"))
            record = {
                "lat": float(lat_i),
                "lon": float(lon_i),
                "i": int(i),
                "j": int(j),
                "data": data
            }
            records.append(record)
        with open(filename, "w") as f:
            json.dump(records, f, indent=2)

    def variables(self, t):
        with netCDF4.Dataset(self.path) as ds:
            # velocities of first row and column are not used
            # TODO: this does not quite match up....
            u1 = np.squeeze(ds.variables['velu'][t][self.s][1:, 1:])
            v1 = np.squeeze(ds.variables['velv'][t][self.s][1:, 1:])
            sep = np.squeeze(ds.variables['sep'][t][self.s][1:, 1:])
        return dict(
            u1=u1,
            v1=v1,
            sep=sep
        )

    def timeseries(self, i, j):
        with netCDF4.Dataset(self.path) as ds:
            u1 = np.squeeze(ds.variables['velu'][:, i, j])
            v1 = np.squeeze(ds.variables['velv'][:, i, j])
            s1 = np.squeeze(ds.variables['sep'][:, i, j])

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

    def validate(self):
        """validate a file"""
        valid = True
        with netCDF4.Dataset(self.path) as ds:
            variables = ds.variables.keys()
        for var in ("velu", "velv", "lat", "lon"):
            if var not in variables:
                logger.warn(
                    "%s not found in variables of file %s",
                    var,
                    self.path
                )
                valid = False
                return valid
        return valid

    def meta(self):
        metadata = super().meta()
        metadata['metadata']['format'] = str(self)
        llur_wgs84 = self.canvas['bbox_wgs84']
        extent = dict(
            # lat, lon
            sw=[llur_wgs84[1], llur_wgs84[0]],
            ne=[llur_wgs84[4], llur_wgs84[3],],
            time=(self.grid['time'][0], self.grid['time'][-1])

        )
        # make sure we have an extent
        metadata['extent'] = metadata.get('extent', {})
        metadata['extent'].update(extent)

        path = pathlib.Path(self.path)
        src = path.with_suffix('.mp4').relative_to(path.parent)
        uv = {
            "src": str(src),
            "type": "video/mp4",
            "tag": "video",
            "width": self.canvas['nx'],
            "height": self.canvas['ny']
        }
        metadata['uv'] = uv
        metadata['uniforms'] = {
            "flipv": True
        }
        return metadata
