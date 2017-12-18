import bisect
import logging

import numpy as np
import scipy.interpolate
import pandas as pd
import tqdm
import geojson

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
    bins = row["bins"]

    # if we don't have any volume table
    if bins is None:
        return None
    volume_table = row["volume_table"]
    cum_volume_table = row["cum_volume_table"]

    dem_i = dem['band'][row['slice']]

    # this part is once volume is known
    vol_i = row['vol1']

    fill_idx = bisect.bisect(cum_volume_table, vol_i)
    remaining_volume = vol_i - cum_volume_table[fill_idx - 1]
    pixel_area = dem['dxp'] * dem['dyp']
    face_area = np.prod(dem_i.shape) * pixel_area

    if fill_idx >= len(cum_volume_table) - 1:
        remaining = (vol_i - cum_volume_table[-1]) / face_area
        target_level = bins[-1] + remaining
    else:
        remaining_volume_fraction = remaining_volume / volume_table[fill_idx]
        target_level = bins[fill_idx] + remaining_volume_fraction * (bins[fill_idx + 1] - bins[fill_idx])
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
    face_centers = grid['face_centers']
    L = scipy.interpolate.LinearNDInterpolator(face_centers, values)
    return L


def build_tables(grid, dem):
    """compute volume tables per cell"""

    # compute cache of histograms per cell
    faces = grid['face_coordinates']
    rows = []
    # TODO: run this in parallel (using concurrent futures)
    for id_, face in tqdm.tqdm(enumerate(faces), total=faces.shape[0], desc='table rows'):
        affine = dem['affine']
        face_px = dem['world2px'](face)
        face_px2slice = np.s_[
            face_px[:, 1].min():face_px[:, 1].max(),
            face_px[:, 0].min():face_px[:, 0].max()
        ]
        dem_i = dem['band'][face_px2slice]
        if dem_i.mask.any():
            n, bins = None, None
            volume_table = None
            cum_volume_table = None
        else:
            n, bins = np.histogram(dem_i, bins=20)
            n_cum = np.cumsum(n)
            volume_table = np.abs(affine.a * affine.e) * n_cum * np.diff(bins)
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
            volume_table=volume_table,
            cum_volume_table=cum_volume_table,
            n=n,
            extent=extent,
            bins=bins
        )
        rows.append(record)


    tables = pd.DataFrame.from_records(rows).set_index('id')
    return tables


def compute_features(dem, tables, data, method='waterdepth'):
    """compute subgrid waterdepth band"""

    # register pandas progress
    tqdm.tqdm(desc="panda is out for lunch!").pandas()

    faces = list(tables.index)

    tables['vol1'] = data['vol1']
    tables['s1'] = data['s1']
    tables['waterdepth'] = data['waterdepth']

    results = []
    # fill the in memory band
    for face_idx in tqdm.tqdm(faces):
        row = tables.loc[face_idx]
        result = subgrid_compute(row, dem=dem, method=method)
        results.append(result)
    tables['subgrid_' + method] = results

    def row2feature(row):
        """convert row 2 features"""
        coordinates = row['face'].mean(axis=0)
        feature = geojson.Feature(
            geometry=geojson.Point(
                coordinates=tuple(coordinates)
            ),
            id=int(row.name),
            properties={
                "s1": float(row.s1),
                "subgrid_" + method: float(row['subgrid_' + method]),
                "vol1": float(row.vol1),
                "waterdepth": float(row.waterdepth)
            }
        )
        return feature
    features = list(
        tables.progress_apply(row2feature, axis=1)
    )
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
    for face_idx in tqdm.tqdm(faces):
        row = tables.loc[face_idx]
        result = subgrid_compute(row, dem=dem, method=method)
        if result is None:
            excluded.append(face_idx)
            continue
        band[row['slice']] = result
    logger.info("skipped %s cells (%s)", len(excluded), excluded)
    return band


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
