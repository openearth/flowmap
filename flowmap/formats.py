# supported file formats
import logging
import functools

import netCDF4
import scipy.interpolate
import numpy as np
import mako.template
import osgeo.osr
import skimage.draw
import skimage.morphology
import matplotlib.pyplot as plt
import matplotlib.colors

matplotlib.use('Agg')

logger = logging.getLogger(__name__)

# what to export on import from *, also used to get a list of available formats
__all__ = ["Matroos", "get_format"]

dump_tmpl = """
% for var in grid:
${var}
 - shape: ${grid[var].shape}
 - type:  ${grid[var].dtype}
 - min:   ${grid[var].min()}
 - max:   ${grid[var].max()}
% endfor

% for var in canvas:
- ${var}: ${canvas[var]}
% endfor


"""


def transform(x, y, transformation):
    """transform coordinates, for n-d coordinates with masks"""
    if len(x.shape) <= 2:
        # curvilinear
        mask = np.logical_or(x.mask, y.mask)
    else:
        # vertices
        mask = np.logical_or(x.mask, y.mask).any(axis=-1)

    xy = np.stack([x[~mask], y[~mask]], axis=-1)
    # shape of the non masked contours
    old_shape = x[~mask].shape
    new_shape = np.prod(xy.shape[:-1]), +  xy.shape[-1],
    xy = xy.reshape(new_shape)
    x_t, y_t, _ = np.array(transformation.TransformPoints(xy)).T
    x_t = x_t.reshape(old_shape)
    y_t = y_t.reshape(old_shape)
    X_t = np.ma.masked_all_like(x)
    Y_t = np.ma.masked_all_like(y)
    X_t[~mask] = x_t
    Y_t[~mask] = y_t
    return X_t, Y_t


def get_format(dataset, **kwargs):
    """get the format for the dataset"""
    for format in [Matroos, Delft3DMatlab]:
        ds = format(dataset, **kwargs)
        if ds.validate():
            logger.info("Found valid format %s for %s", format, dataset)
            return format
    return Matroos


class NetCDF(object):
    def __init__(self, path, src_epsg=4326, dst_epsg=28992):
        self.path = path
        # source and destination epsg code
        self.src_epsg = src_epsg
        self.dst_epsg = dst_epsg

    @property
    def srs(self):
        # let's define the systems
        src_srs = osgeo.osr.SpatialReference()
        src_srs.ImportFromEPSG(self.src_epsg)

        # google mercator
        web_srs = osgeo.osr.SpatialReference()
        web_srs.ImportFromEPSG(3857)

        # Lat,Lon
        wgs84 = osgeo.osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)

        # local UTM
        utm = osgeo.osr.SpatialReference()
        utm.ImportFromEPSG(self.dst_epsg)

        # and the translations between them
        src2wgs84 = osgeo.osr.CoordinateTransformation(src_srs, wgs84)
        web2wgs84 = osgeo.osr.CoordinateTransformation(web_srs, wgs84)
        utm2wgs84 = osgeo.osr.CoordinateTransformation(utm, wgs84)
        wgs842utm = osgeo.osr.CoordinateTransformation(wgs84, utm)
        wgs842web = osgeo.osr.CoordinateTransformation(wgs84, web_srs)
        utm2web = osgeo.osr.CoordinateTransformation(utm, web_srs)
        src2utm = osgeo.osr.CoordinateTransformation(src_srs, utm)
        src2web = osgeo.osr.CoordinateTransformation(src_srs, web_srs)

        return dict(
            src2wgs84=src2wgs84,
            web2wgs84=web2wgs84,
            utm2wgs84=utm2wgs84,
            wgs842utm=wgs842utm,
            wgs842web=wgs842web,
            utm2web=utm2web,
            src2utm=src2utm,
            src2web=src2web
        )

    def dump(self):
        tmpl = mako.template.Template(dump_tmpl)
        text = tmpl.render(grid=self.grid, canvas=self.canvas)
        return text


class Delft3DMatlab(NetCDF):
    """NetCDF converted with vs_nefis2nc"""

    @property
    @functools.lru_cache()
    def grid(self):
        """generate grid variables"""
        src2wgs84 = self.srs['src2wgs84']
        src2utm = self.srs['src2utm']
        src2web = self.srs['src2web']
        utm2web = self.srs['utm2web']

        grid = {}
        with netCDF4.Dataset(self.path) as ds:
            times = netCDF4.num2date(ds.variables['time'][:],
                                     ds.variables['time'].units)

            grid['times'] = times
            x, y = ds.variables['x'][:], ds.variables['y'][:]
            mask = np.logical_or(x.mask, y.mask)

            grid['x'] = x
            grid['y'] = y

            Lon, Lat = transform(x, y, src2wgs84)
            X_web, Y_web = transform(x, y, src2web)
            X_utm, Y_utm = transform(x, y, src2utm)

            # return
            grid["Lon"] = Lon
            grid["Lat"] = Lat
            grid["X_web"] = X_web
            grid["Y_web"] = Y_web
            grid["X_utm"] = X_utm
            grid["Y_utm"] = Y_utm

            # same thing for contours
            xcc = ds.variables['grid_x'][:]
            ycc = ds.variables['grid_y'][:]

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

            X_utm_c_a, X_utm_c_b = transform(xcc, ycc, src2utm)

            grid['X_distort'], grid['Y_distort'] = self.compute_distortion(X_utm, Y_utm, utm2web)

            # initial values (used to determine shapes and stuff, maybe remove if not used)
            grid["s1_0"] = np.squeeze(ds.variables['waterlevel'][0])
            grid["u1_0"] = np.squeeze(ds.variables['velocity_x'][0])
            grid["v1_0"] = np.squeeze(ds.variables['velocity_y'][0])

        return grid

    def animate(self):
        """generate an animation (set of png's)"""

        F, X, Y = self.canvas['F'], self.canvas['X'], self.canvas['Y']

        X_web = self.grid['X_web']
        Y_web = self.grid['Y_web']
        coordinate_mask = np.logical_or(X_web.mask, Y_web.mask)

        vmin, vmax = -0.2, 0.2
        N = matplotlib.colors.Normalize(vmin, vmax, clip=True)
        for i, t in enumerate(self.grid['time']):
            # get data for t = t
            data = self.variables(i)
            # compute uv and mask with coordinate mask
            uv = np.c_[data['u1'][~coordinate_mask], data['v1'][~coordinate_mask].ravel()]
            F.values = uv.astype(F.values.dtype)
            UV = F(X, Y)
            RG = N(UV)
            R, G = RG[..., 0], RG[..., 1]

            # cells without a velocity
            value_mask = np.logical_and(UV[..., 0] == 0.0, UV[..., 1] == 0.0)
            # masked cells
            B = np.zeros_like(R) + np.logical_and(self.canvas['is_grid'], ~value_mask)
            RGB = np.dstack([R, G, B])
            # store in filename
            # TODO: generate with ffmpeg
            plt.imsave('test_%04d.png' % (i,), RGB)

    def variables(self, t):
        with netCDF4.Dataset(self.path) as ds:
            u1 = np.squeeze(ds.variables['velocity_x'][t])
            v1 = np.squeeze(ds.variables['velocity_y'][t])
        return dict(
            u1=u1,
            v1=v1
        )

    @property
    def canvas(self):
        """determine the rendering canvas and compute coordinates"""
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


class Matroos(NetCDF):
    """FEWS Matroos format"""
    @property
    def grid(self):
        """generate global variables"""
        with netCDF4.Dataset(self.path) as ds:
            times = netCDF4.num2date(ds.variables['time'][:],
                                     ds.variables['time'].units)
            analysis_time = netCDF4.num2date(ds.variables['analysis_time'][:],
                                             ds.variables['analysis_time'].units)

            lat = ds.variables['lat'][:]
            lon = ds.variables['lon'][:]

            # initial values (used to determine shapes and stuff, maybe remove if not used)
            sep_0 = ds.variables['sep'][0]
            u1_0 = ds.variables['velu'][0]
            v1_0 = ds.variables['velv'][0]

        variables = dict(
            times=times,
            analysis_time=analysis_time,
            u1_0=u1_0,
            v1_0=v1_0,
            sep_0=sep_0,
            lat=lat,
            lon=lon
        )
        return variables

    @property
    def canvas(self):
        return {}

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
