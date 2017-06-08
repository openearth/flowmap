# supported file formats
import logging
import functools

import numpy as np
import matplotlib.colors
import tqdm

matplotlib.use('Agg')

logger = logging.getLogger(__name__)

# what to export on import from *, also used to get a list of available formats


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
    # make sure you don't use a masked array (very slooow)
    xy = xy.reshape(new_shape).filled()
    xy_t = transformation.TransformPoints(xy)
    x_t, y_t, _ = np.array(xy_t).T
    x_t = x_t.reshape(old_shape)
    y_t = y_t.reshape(old_shape)
    X_t = np.ma.masked_all_like(x)
    Y_t = np.ma.masked_all_like(y)
    X_t[~mask] = x_t
    Y_t[~mask] = y_t
    return X_t, Y_t


def get_format(dataset, **kwargs):
    """get the format for the dataset"""

    from .matroos import Matroos
    from .delft3d import Delft3DMatlab

    for format in [Matroos, Delft3DMatlab]:
        ds = format(dataset, **kwargs)
        if ds.validate():
            logger.info("Found valid format %s for %s", format, dataset)
            return format
    return Matroos


@functools.lru_cache()
def build_fill_rules():
    """define fill rules for a 5x5 matrix"""
    m = 2
    rules = []
    for (i, j) in [
            # + shape
            (m, m),
            (m, m-1),
            (m+1, m),
            (m, m+1),
            (m-1, m),
            (m-1, m-1),
            (m+1, m-1),
            (m-1, m+1),
            (m+1, m+1)
    ]:
        rule = {
            "ij": [i, j],
            "interpolates": [
                [[i-1, j], [i+1, j]],
                [[i, j-1], [i, j+1]]
            ],
            "extrapolate": []
        }
        if i >= 2:
            rule["extrapolate"].append([[i-1, j], [i-2, j]])
        if i < 3:
            rule["extrapolate"].append([[i+1, j], [i+2, j]])
        if j >= 2:
            rule["extrapolate"].append([[i, j-1], [i, j-2]])
        if j < 3:
            rule["extrapolate"].append([[i, j+1], [i, j+2]])
        rules.append(rule)
    return rules


def window(i, j, lat, lon):
    """create a padded window of 5x5 around point i,j for the arrays lat and lon """
    m, n = lat.shape
    i_min = max(i-2, 0)
    i_max = min(i+3, m)
    j_min = max(j-2, 0)
    j_max = min(j+3, n)
    s = np.s_[i_min:i_max, j_min:j_max]
    views = []
    for arr in [lat, lon]:
        view = arr[s].copy()
        if view.shape == (5, 5):
            views.append(view)
        else:
            i_0 = max(2 - i, 0)
            j_0 = max(2 - j, 0)
            padded_view = np.ma.masked_all((5, 5), dtype=view.dtype)
            padded_view[i_0:i_0+view.shape[0], j_0:j_0+view.shape[1]] = view
            views.append(padded_view)
    # should contain lat_i, lon_i
    return tuple(views)


def ij2verts(i, j, lat, lon, rules):
    """compute cell boundaries for points """
    lat_i, lon_i = window(i, j, lat, lon)
    for obj in rules:
        i, j = obj['ij']
        if not lat_i.mask[i, j]:
            # should not be masked
            continue
        # filling with obj
        weight = 0
        lat_ij = 0
        lon_ij = 0
        for (i_a, j_a), (i_b, j_b) in obj["interpolates"]:
            # can we interpolate
            if (not lat_i.mask[i_a, j_a]) and (not lat_i.mask[i_b, j_b]):
                lat_ij += lat_i[i_a, j_a]
                lat_ij += lat_i[i_b, j_b]
                lon_ij += lon_i[i_a, j_a]
                lon_ij += lon_i[i_b, j_b]
                weight += 2.0
        if weight > 0:
            # interpolate
            lat_i[i, j] = lat_ij/weight
            lon_i[i, j] = lon_ij/weight
            # this point is filled, on with the next one
            continue
        # if we can't interpolate let's try to extrapolate
        weight = 0
        lat_ij = 0
        lon_ij = 0
        for (i_a, j_a), (i_b, j_b) in obj["extrapolate"]:
            # can we extrapolate
            if (not lat_i.mask[i_a, j_a]) and (not lat_i.mask[i_b, j_b]):
                lon_ij += lon_i[i_a, j_a] + (lon_i[i_a, j_a] - lon_i[i_b, j_b])
                lat_ij += lat_i[i_a, j_a] + (lat_i[i_a, j_a] - lat_i[i_b, j_b])
                weight += 1
        if weight > 0:
            lat_i[i, j] = lat_ij/weight
            lon_i[i, j] = lon_ij/weight
            # this point is filled, on with the next one
            continue
    m = 2

    ll_idx = np.array([[m, m], [m, m-1], [m+1, m-1], [m+1, m]])
    lr_idx = np.array([[m, m], [m+1, m], [m+1, m+1], [m, m+1]])
    ur_idx = np.array([[m, m], [m, m+1], [m-1, m+1], [m-1, m]])
    ul_idx = np.array([[m, m], [m-1, m], [m-1, m-1], [m, m-1]])

    vertices = []
    for idx in (ll_idx, lr_idx, ur_idx, ul_idx):
        s = idx[:, 0], idx[:, 1]

        vertex = (lon_i[s], lat_i[s])
        if np.ma.is_masked(vertex[0]) or np.ma.is_masked(vertex[1]):
            continue
        vertices.append(vertex)
    verts = np.array(vertices)
    return verts


def points2contours(lat, lon):
    """convert curvilinear grid defined by lat and lon to cell contours"""
    # define
    rules = build_fill_rules()
    vertices = {}
    for i, j in tqdm.tqdm(list(np.ndindex(lat.shape)), desc="contours"):
        verts = ij2verts(i, j, lat, lon, rules)
        if len(verts) and verts.shape == (4, 2, 4):
            vertices[(i, j)] = verts.mean(axis=-1)
    return vertices


def contours2vertices(lat, lon):
    """convert curvilinear contour grid defined by lat and lon to cell contours"""
    # lookup cell borders ccw
    ul = np.ma.vstack([lon[:-1, :-1].ravel(), lat[:-1, :-1].ravel()]).T[:, np.newaxis, :]
    ur = np.ma.vstack([lon[:-1, 1:].ravel(), lat[:-1, 1:].ravel()]).T[:, np.newaxis, :]
    lr = np.ma.vstack([lon[1:, 1:].ravel(), lat[1:, 1:].ravel()]).T[:, np.newaxis, :]
    ll = np.ma.vstack([lon[1:, :-1].ravel(), lat[1:, :-1].ravel()]).T[:, np.newaxis, :]

    vertices = np.ma.hstack([ul, ur, lr, ll])
    # is any of the coords missing
    mask = vertices.mask.any(axis=(1, 2))
    # reset mask
    vertices.mask[:] = False
    vertices.mask[mask] = True
    return vertices
