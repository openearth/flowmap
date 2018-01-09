import numpy as np

import networkx
from tvtk.api import tvtk
from tvtk.common import configure_input, configure_source_data


def concave_hull(polydata):
    """extract the concave hull of a vtk polydata object"""
    # setup the extraction pipeline
    extract = tvtk.FeatureEdges()
    # we're only interested in the boundaries
    extract.feature_edges = False
    extract.boundary_edges = True
    extract.manifold_edges = False
    extract.non_manifold_edges = False
    configure_input(extract, polydata)
    # compute edges
    extract.update()
    # extract the points
    points = extract.output.points.to_array()

    # slice every 2nd and 3rd point
    line_ids = np.c_[
        extract.output.lines.to_array()[1::3],
        extract.output.lines.to_array()[2::3]
    ]
    # 1st points should all be 2
    assert (extract.output.lines.to_array()[0::3] == 2).all(), "expected only lines"

    # construct a directed graph
    D = networkx.DiGraph(data=line_ids.tolist())
    # find the first cycle
    first_cycle = networkx.find_cycle(D)
    # create the index to lookup points, including last point
    cycle = [x[0] for x in first_cycle] + [first_cycle[0][0]]
    # fill in index and return first 2 coordinates
    return points[cycle, :2]
