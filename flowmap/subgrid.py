import bisect
import logging

import numpy as np
import pandas as pd
import tqdm
import geojson
import netCDF4
import rasterio
import osgeo.osr
import osgeo.gdal

logger = logging.getLogger(__name__)

# use the same nodata everywhere
NODATA = -9999


class MetaArray(np.ndarray):
    """Array with metadata."""
    def __new__(cls, array, dtype=None, order=None, **kwargs):
        obj = np.asarray(array, dtype=dtype, order=order).view(cls)
        obj.metadata = kwargs
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.metadata = getattr(obj, 'metadata', None)


def data_for_idx(face_idx, dem, grid, data):
    """get data for cell with face_idx face_idx"""
    face = grid['face_coordinates'][face_idx]
    affine = dem['affine']
    idx = (face - (affine.xoff, affine.yoff)) / (affine.a, affine.e)
    i_min, i_max = int(idx[:, 0].min()), int(idx[:, 0].max())
    j_min, j_max = int(idx[:, 1].min()), int(idx[:, 1].max())
    dem_i = dem['band'][j_min:j_max, i_min:i_max]
    vol_i = data['vol1'][face_idx]
    data = dict(
        face=face,
        dem=dem_i,
        vol=vol_i
    )
    return data


def compute_waterlevel_per_cell(row, dem):
    """get the subgrid waterdepth image"""
    # tables is a dataframe
    bin_edges = row["bin_edges"]

    # if we don't have any volume table
    if bin_edges is None:
        return None
    volume_table = row["volume_table"]
    cum_volume_table = row["cum_volume_table"]

    # this part is once volume is known
    vol_i = row['vol1']

    fill_idx = bisect.bisect(cum_volume_table, vol_i)
    if fill_idx > 0:
        remaining_volume = vol_i - cum_volume_table[fill_idx - 1]
    else:
        remaining_volume = vol_i - cum_volume_table[0]
    face_area = row['face_area']

    # are we outside the volume table
    if fill_idx >= len(cum_volume_table) - 1:
        remaining = (vol_i - cum_volume_table[-1]) / face_area
        target_level = bin_edges[-1] + remaining
    else:
        # we're in the volume table
        remaining_volume_fraction = remaining_volume / volume_table[fill_idx]
        target_level = bin_edges[fill_idx] + remaining_volume_fraction * (bin_edges[fill_idx + 1] - bin_edges[fill_idx])
    result = float(target_level)
    return result


def build_tables(ugrid, dem, id_grid, valid_range=None):
    """compute volume tables per cell"""

    if (valid_range is not None) and (None not in valid_range):
        logger.info('filtering by valid-range %s', valid_range)
        invalid_mask = np.logical_or(
            dem['band'] < valid_range[0],
            dem['band'] > valid_range[1]
        )
    else:
        invalid_mask = np.zeros_like(dem['band'], dtype=bool)

    # compute cache of histograms per cell
    faces = ugrid['face_coordinates']
    rows = []
    # TODO: run this in parallel (using concurrent futures)
    for id_, face in tqdm.tqdm(enumerate(faces), total=faces.shape[0], desc='table rows'):
        # # Use this for faster debugging of triangles
        # if id_ < 4700:
        #     continue
        # if id_ > 4800:
        #     break

        affine = dem['affine']
        # remove masked coordinates
        face = face[~face.mask.any(axis=1), :]
        face_px = dem['world2px'](face)
        face_px2slice = np.s_[
            face_px[:, 1].min():face_px[:, 1].max(),
            face_px[:, 0].min():face_px[:, 0].max()
        ]
        face_slice = [
            face_px2slice[0].start,
            face_px2slice[0].stop,
            face_px2slice[1].start,
            face_px2slice[1].stop
        ]
        dem_i = dem['band'][face_px2slice]
        ids_i_mask = id_grid[face_px2slice] != id_
        # we have three conditions to exclude a cell
        masks = [
            # dem is missing
            dem_i.mask,
            # not our cell
            ids_i_mask,
        ]
        # cell not set
        if hasattr(ids_i_mask, 'mask') and ids_i_mask.mask.any():
            masks.append(ids_i_mask.mask)
        # invalid values
        if valid_range is not None:
            masks.append(
                invalid_mask[face_px2slice]
            )

        mask = np.logical_or.reduce(masks)
        if dem_i.mask.any() or mask.all():
            n_per_bin, bin_edges = None, None
            volume_table = None
            cum_volume_table = None
            face_area = None
        else:
            n_per_bin, bin_edges = np.histogram(dem_i[~mask], bins=20)
            # should this be equal to non masked cells in dem_i?
            n_cum = np.cumsum(n_per_bin)
            volume_table = np.abs(affine.a * affine.e) * n_cum * np.diff(bin_edges)
            face_area = np.abs(affine.a * affine.e) * n_cum
            cum_volume_table = np.cumsum(volume_table)
        extent = [
            face[:, 0].min(),
            face[:, 0].max(),
            face[:, 1].min(),
            face[:, 1].max()
        ]
        record = dict(
            id=id_,
            slice=face_slice,
            face=face,
            face_area=face_area,
            volume_table=volume_table,
            cum_volume_table=cum_volume_table,
            n_per_bin=n_per_bin,
            extent=extent,
            bin_edges=bin_edges
        )
        rows.append(record)

    tables = pd.DataFrame.from_records(rows).set_index('id')
    return tables


def compute_waterlevels(grid, dem, tables, data):
    """compute subgrid waterdepth band"""

    # register pandas progress
    tqdm.tqdm(desc="computing features").pandas()

    # list of face indices
    face_ids = np.arange(tables['volume_table'].shape[0])

    tables['vol1'] = data['vol1']
    tables['s1'] = data['s1']
    tables['waterdepth'] = data['waterdepth']

    results = []
    # fill the in memory band
    for face_id in tqdm.tqdm(face_ids, desc='subgrid compute'):
        row = {}
        for key, var in tables.items():
            if key == 'metadata':
                continue

            row[key] = var[face_id]
        result = compute_waterlevel_per_cell(row, dem=dem)
        results.append(result)
    tables['subgrid_waterlevel'] = results

    features = []
    centroids = grid['face_centroids']
    for face_id in tqdm.tqdm(face_ids, desc='exporting features'):
        """convert row 2 features"""
        centroid = centroids[face_id]
        feature = geojson.Feature(
            geometry=geojson.Point(
                coordinates=tuple(centroid)
            ),
            id=face_id,
            properties={
                "s1": data['s1'][face_id],
                "subgrid_waterlevel": tables['subgrid_waterlevel'][face_id],
                "vol1": tables['vol1'][face_id],
                "waterdepth": tables['waterdepth'][face_id]
            }
        )
        features.append(feature)
    collection = geojson.FeatureCollection(features=features)
    return collection


def interpolate_waterlevels(waterlevel_file, interpolated_waterlevel_file, dem, epsg):
    """Compute the interpolated waterlevels in the same format as dem. Reults are stored in interpolated_waterlevel_file."""
    # compute bounds from dem
    ulx, xres, xskew, uly, yskew, yres = dem['transform']
    lrx = ulx + (dem['width'] * xres)
    lry = uly + (dem['height'] * yres)
    # define CRS
    srs = osgeo.osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    # create options object
    options = osgeo.gdal.GridOptions(
        zfield='subgrid_waterlevel',
        outputSRS=srs,
        outputType=osgeo.gdal.GDT_Float32,
        noData=NODATA,
        width=dem['width'],
        height=dem['height'],
        outputBounds=(ulx, uly, lrx, lry),
        algorithm='invdistnn:power=3.0:max_points=4:radius=8:nodata={}'.format(NODATA)
    )

    logger.info('interpolating')
    # note dst, src order
    result = osgeo.gdal.Grid(
        str(interpolated_waterlevel_file),
        str(waterlevel_file),
        options=options
    )
    result.FlushCache()
    return


def create_export(filename, n_cells, n_bins, attributes=None):
    """create an export file for subgrid tables"""

    dimensions = {
        "cells": n_cells,
        "bins": n_bins,
        "bin_edges": n_bins + 1,
        "two_times_two": 4
    }
    variables = [
        {
            "name": "bin_edges",
            "dimensions": ("cells", "bin_edges"),
            "long_name": "bin edges of topography histogram",
            "type": "double"
        },
        {
            "name": "cum_volume_table",
            "dimensions": ("cells", "bins"),
            "long_name": "cumulative volume table",
            "type": "double"
        },
        {
            "name": "volume_table",
            "dimensions": ("cells", "bins"),
            "long_name": "volume table",
            "type": "double"
        },
        {
            "name": "face_area",
            "dimensions": ("cells", ),
            "long_name": "face area",
            "type": "double"
        },
        {
            "name": "extent",
            "dimensions": ("cells", "two_times_two"),
            "long_name": "extent (left, right, lower, upper)",
            "type": "double"
        },
        {
            "name": "n_per_bin",
            "dimensions": ("cells", "bins"),
            "long_name": "topography histogram",
            "type": "int"
        },
        {
            "name": "slice",
            "dimensions": ("cells", "two_times_two"),
            "long_name": "slice (row start, stop, colum start stop)",
            "type": "int"
        }
    ]

    with netCDF4.Dataset(filename, 'w') as ds:
        for name, size in dimensions.items():
            ds.createDimension(name, size)
        for var in variables:
            ncvar = ds.createVariable(
                var['name'],
                datatype=var['type'],
                dimensions=var['dimensions']
            )
            ncvar.setncattr('long_name', var['long_name'])
        if attributes:
            # store attribute
            for key, val in attributes.items():
                try:
                    ds.setncattr(key, val)
                except TypeError:
                    logger.exception('could not set {} to {}'.format(key, val))

def export_tables(filename, tables):
    """store tables in netcdf file, create file with create_export"""
    with netCDF4.Dataset(filename, 'r+') as ds:
        for (i, row) in tqdm.tqdm(
            tables.reset_index().iterrows(),
            total=len(tables),
            desc='exporting'
        ):
            for var in [
                'bin_edges', 'cum_volume_table', 'volume_table',
                'extent', 'n_per_bin', 'face_area', 'slice'
            ]:
                val = row[var]
                # skip none
                if val is None:
                    continue

                ds.variables[var][i] = val


def import_tables(filename):
    """import tables from netcdf table dump"""
    with netCDF4.Dataset(filename) as ds:
        # lookup metadata
        metadata = {
            key: getattr(ds, key)
            for key
            in ds.ncattrs()
        }
        vars = {}
        for var in [
            'bin_edges', 'cum_volume_table',
            'volume_table', 'extent', 'n_per_bin',
            'slice', 'face_area'
        ]:
            arr = ds.variables[var][:]
            vars[var] = arr
    # fast approach, use just numpy array
    tables = vars
    tables['metadata'] = metadata
    return tables


def build_id_grid(polys, dem):
    """create a map in the same shape as dem, with face number for each pixel"""
    # generate geojson polygons
    # convert to raster with id as property
    rasterized = rasterio.features.rasterize(
        ((poly, i) for (i, poly) in enumerate(polys)),
        out_shape=dem['band'].shape,
        transform=dem['affine'],
        fill=NODATA
    )
    return rasterized


def export_grid(filename, grid, affine, width, height, epsg):
    """export the id grid to filename"""
    options = dict(
        dtype=str(grid.dtype),
        nodata=NODATA,
        count=1,
        compress='lzw',
        tiled=True,
        blockxsize=256,
        blockysize=256,
        driver='GTiff',
        affine=affine,
        width=width,
        height=height,
        crs=rasterio.crs.CRS({'init': 'epsg:%d' % (epsg, )})
    )

    # fill missings if used (not quite sure if needed)
    # TODO: check if we can remove this.
    if hasattr(grid, 'filled'):
        grid = grid.filled(NODATA)

    with rasterio.open(str(filename), 'w', **options) as out:
        out.write(grid, indexes=1)


def import_id_grid(id_grid_name):
    with rasterio.open(str(id_grid_name)) as src:
        # read band 0 (1-based)
        band = src.read(1, masked=True)
        return band
