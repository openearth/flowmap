import bisect
import logging

import numpy as np
import scipy.interpolate
import pandas as pd
import tqdm
import geojson
import netCDF4

logger = logging.getLogger(__name__)


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


def subgrid_compute(row, dem, method="waterlevel"):
    """get the subgrid waterdepth image"""
    # tables is a dataframe
    bin_edges = row["bin_edges"]

    # if we don't have any volume table
    if bin_edges is None:
        return None
    volume_table = row["volume_table"]
    cum_volume_table = row["cum_volume_table"]
    # imported data (creating slice objects is a bit slow)
    dem_i = dem['band'][
        row['slice'][0]:row['slice'][1],
        row['slice'][2]:row['slice'][3]
    ]

    # this part is once volume is known
    vol_i = row['vol1']

    fill_idx = bisect.bisect(cum_volume_table, vol_i)
    if fill_idx > 0:
        remaining_volume = vol_i - cum_volume_table[fill_idx - 1]
    else:
        remaining_volume = vol_i - cum_volume_table[0]
    pixel_area = dem['dxp'] * dem['dyp']

    # face_area = np.prod(dem_i.shape) * pixel_area
    # TODO: check if this works
    face_area = row['face_area']

    # are we outside the volume table
    if fill_idx >= len(cum_volume_table) - 1:
        remaining = (vol_i - cum_volume_table[-1]) / face_area
        target_level = bin_edges[-1] + remaining
    else:
        # we're in the volume table
        remaining_volume_fraction = remaining_volume / volume_table[fill_idx]
        target_level = bin_edges[fill_idx] + remaining_volume_fraction * (bin_edges[fill_idx + 1] - bin_edges[fill_idx])
    if method == 'waterlevel':
        result = float(target_level)
    elif method == 'waterdepth':
        # first cell that is not completely filled
        waterdepth_i = np.zeros_like(dem_i)
        idx = dem_i < target_level
        waterdepth_i[idx] = (target_level - dem_i[idx])
        result = waterdepth_i
    return result


def build_interpolate(grid, values):
    """create an interpolation function"""
    # assert a pyugrid
    face_centroids = grid['face_centroids']
    L = scipy.interpolate.LinearNDInterpolator(face_centroids, values)
    return L


def build_tables(grid, dem, id_grid, valid_range):
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
    faces = grid['face_coordinates']
    rows = []
    # TODO: run this in parallel (using concurrent futures)
    for id_, face in tqdm.tqdm(enumerate(faces), total=faces.shape[0], desc='table rows'):
        # Use this for faster debugging of triangles
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
        dem_i = dem['band'][face_px2slice]
        ids_i_mask  = id_grid[face_px2slice] != id_
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
        # TOOD: als mask using id grid here....
        if dem_i.mask.any():
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
            slice=face_px2slice,
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


def compute_features(grid, dem, tables, data, method='waterdepth'):
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
            row[key] = var[face_id]

        result = subgrid_compute(row, dem=dem, method=method)
        results.append(result)
    tables['subgrid_' + method] = results

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
                "subgrid_" + method: tables['subgrid_' + method][face_id],
                "vol1": tables['vol1'][face_id],
                "waterdepth": tables['waterdepth'][face_id]
            }
        )
        features.append(feature)
    collection = geojson.FeatureCollection(features=features)
    return collection


def compute_band(grid, dem, tables, data, method='waterdepth'):
    """compute subgrid waterdepth band"""
    excluded = []
    faces = list(tables.index)

    # create a masked array, always return floats (band is sometimes int)
    band = np.ma.masked_all(dem['band'].shape, dtype='float')

    tables['vol1'] = data['vol1']

    # fill the in memory band
    for face_idx in tqdm.tqdm(faces, desc='computing band'):
        row = tables.loc[face_idx]
        result = subgrid_compute(row, dem=dem, method=method)
        if result is None:
            excluded.append(face_idx)
            continue
        band[row['slice']] = result
    logger.info("skipped %s cells (%s)", len(excluded), excluded)
    return band


def create_export(filename, n_cells, n_bins):
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


def export_tables(filename, tables):
    """store tables in netcdf file, create file with create_export"""
    with netCDF4.Dataset(filename, 'r+') as ds:
        for i, row in tqdm.tqdm(
                tables.reset_index().iterrows(),
                total=len(tables),
                desc='exporting'
        ):
            for var in [
                'bin_edges', 'cum_volume_table', 'volume_table',
                'extent', 'n_per_bin', 'face_area'
            ]:
                val = row[var]
                # skip none
                if val is None:
                    continue

                ds.variables[var][i] = val
            ds.variables['slice'][i] = [
                row.slice[0].start,
                row.slice[0].stop,
                row.slice[1].start,
                row.slice[1].stop
            ]


def import_tables(filename, arrays=True):
    """import tables from netcdf table dump"""
    with netCDF4.Dataset(filename) as ds:
        vars = {}
        index = np.arange(ds.variables['bin_edges'].shape[0])
        for var in [
                'bin_edges', 'cum_volume_table',
                'volume_table', 'extent', 'n_per_bin',
                'slice', 'face_area'
        ]:
            arr = ds.variables[var][:]
            vars[var] = arr
    # fast approach, use just numpy array
    if arrays:
        tables = vars
    else:
        # this is slow for large number of cells
        for key, arr in vars.items():
            vars[key] = list(arr)
        tables = pd.DataFrame(vars, index=index)
    return tables


def compute_interpolated(L, dem, data, s=None):
    """compute a map of interpolated waterdepth, masked where detailed topography >= interpolated waterlevel, optionally sliced by a tuple (s) of row, column slices"""
    if s is None:
        s = np.s_[:, :]

    # create the pixel grid (assuming no rotation)
    affine = dem['affine']
    assert affine.b == 0 and affine.d == 0, 'rotated dems not implemented'
    y = np.arange(affine.f, affine.f + affine.e * dem['height'], affine.e)
    x = np.arange(affine.c, affine.c + affine.a * dem['width'], affine.a)
    # we need the full grid to get the interpolated values
    X, Y = np.meshgrid(x[s[1]], y[s[0]])
    # fill the interpolation function
    msg = 'Interpolation function should be filled with s1, vol1, and waterdepth'
    assert L.values.shape[1] == 3, msg
    # fill in new values
    L.values = np.c_[data['s1'], data['vol1'], data['waterdepth']]
    # compute interplation
    interpolated = L(X, Y)
    # get the variables
    s1 = interpolated[..., 0]
    waterdepth = interpolated[..., 2]
    vol1 = interpolated[..., 1]
    # lookup band
    dem_band = dem['band'][s]
    # mask interpolated values using dem
    masked_waterdepth = np.ma.masked_array(waterdepth, mask=dem_band >= s1)
    return {
        "masked_waterdepth": masked_waterdepth,
        "s1": s1,
        "vol1": vol1,
        "dem": dem_band
    }
