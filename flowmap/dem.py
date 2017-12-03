"""
Functions related to reading dem files
"""
import rasterio
import rasterio.plot
import numpy as np


def read_dem(dem_filename):
    """read dem"""
    dem = {}
    with rasterio.open(str(dem_filename)) as src:
        # read band 0 (1-based)
        band = src.read(1, masked=True)
        dem['band'] = band
    # pixel sizes
    affine = src.affine
    dem['dxp'] = affine.a
    dem['dyp'] = -affine.e
    dem['affine'] = affine
    dem['width'] = src.width
    dem['height'] = src.height
    def world2px(xy):
        xy_t = (xy - (affine.xoff, affine.yoff)) / (affine.a, affine.e)
        return xy_t

    dem['world2px'] = lambda xy: np.vstack((~dem['affine']) * (xy[:,0], xy[:,1])).T.astype('int')
    dem['px2world'] = lambda xy: np.vstack((dem['affine']) * (xy[:,0], xy[:,1])).T

    # extent for imshow
    dem['extent'] = rasterio.plot.plotting_extent(src)
    return dem
