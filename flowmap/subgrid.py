import bisect
import logging

import numpy as np
import pandas as pd
import tqdm

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


def subgrid_waterdepth(face_idx, dem, grid, data, tables):
    """get the subgrid waterdepth image"""
    # tables is a dataframe
    hist = tables.loc[face_idx]
    bins = hist["bins"]

    # if we don't have any volume table
    if bins is None:
        return None
    volume_table = hist["volume_table"]
    cum_volume_table = hist["cum_volume_table"]

    dem_i = hist['dem']

    # this part is once volume is known
    vol_i = data['vol1'][face_idx]

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

    # first cell that is not completely filled
    waterdepth_i = np.zeros_like(dem_i)
    idx = dem_i < target_level
    waterdepth_i[idx] = (target_level - dem_i[idx])
    return waterdepth_i


def build_tables(grid, dem):
    """compute volume tables per cell"""

    # compute cache of histograms per cell
    idx = range(grid['face_coordinates'].shape[0])
    faces = grid['face_coordinates'][idx]
    tables = {}
    for id_, face in zip(idx, tqdm.tqdm(faces)):
        affine = dem['affine']
        face_px = dem['world2px'](face)
        face_px2slice = np.s_[
            face_px[:, 1].min():face_px[:, 1].max(),
            face_px[:, 0].min():face_px[:, 0].max()
        ]
        dem_i = dem['band'][face_px2slice]
        if not dem_i.mask.all():
            n, bins = np.histogram(dem_i.ravel(), bins=20)
            n_cum = np.cumsum(n)
            volume_table = np.abs(affine.a * affine.e) * n_cum * np.diff(bins)
            cum_volume_table = np.cumsum(volume_table)
        else:
            n, bins = None, None
            volume_table = None,
            cum_volume_table = None
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
            dem=dem_i,
            volume_table=volume_table,
            cum_volume_table=cum_volume_table,
            n=n,
            extent=extent,
            bins=bins
        )
        tables[id_] = record

    tables = pd.DataFrame.from_records(list(tables.values())).set_index('id')
    return tables


def compute_band(grid, dem, tables, data):
    """compute subgrid waterdepth band"""
    excluded = []
    faces = list(tables.index)

    # create a masked array
    band = np.ma.masked_all_like(dem['band'])

    # fill the in memory band
    for face_idx in tqdm.tqdm_notebook(faces):
        row = tables.loc[face_idx]
        waterdepth_i = subgrid_waterdepth(face_idx, dem=dem, grid=grid, data=data, tables=tables)
        if waterdepth_i is None:
            excluded.append(face_idx)
            continue
        band[row['slice']] = waterdepth_i
    logger.info("skipped %s cells (%s)", len(excluded), excluded)
    return band
