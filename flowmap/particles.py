import numpy as np
import pandas as pd
from tvtk.api import tvtk
from tvtk.common import configure_input, configure_source_data, configure_input_data
import vtk
import geojson
import pyproj

def make_tracer_pipeline(polydata, seed):
    # create elements of the pipeline
    tracer = tvtk.StreamTracer()


    # # # You can compute up to 1km per timestep (8km/hr for )

    n = 200  # max number of steps
    l = 1000


    # maximum 1km
    tracer.maximum_propagation = l
    # # # In m
    tracer.integration_step_unit = vtk.vtkStreamTracer.LENGTH_UNIT
    # # # Minimum 5 per step
    tracer.minimum_integration_step = (l/n)
    # # # Maximum 100m per step
    tracer.maximum_integration_step = 10*(l/n)
    # # # Maximum 200 steps
    tracer.maximum_number_of_steps = n
    # # # Maximum error 1cm
    tracer.maximum_error = 1e-2
    # # # We use a path integration. You could argue that you need a
    # # # particle tracking algorithm that matches the numerical grid
    # # # (in our case edge velocities
    # # # and integration over a cell instead of over a line)
    tracer.integrator_type = 'runge_kutta45'

    tracer.debug = True

    cell2point = tvtk.CellDataToPointData()
    # setup the pipeline
    configure_input(cell2point, polydata)
    configure_input(tracer, cell2point)
    configure_source_data(tracer, seed)
    return tracer

def extract_points(tracer):
    points = tracer.output.points.to_array()
    points_data = {
        "point_idx": np.arange(points.shape[0]),
        "x": points[:, 0],
        "y": points[:, 1],
        "z": points[:, 2]
    }

    for i in range(tracer.output.point_data.number_of_arrays):
        name = tracer.output.point_data.get_array_name(i)
        arr = tracer.output.point_data.get_array(i).to_array()
        points_data[name] = list(arr)
    return pd.DataFrame(points_data)

def extract_lines(tracer):
    """convert the streamlines to a data frame"""
    data = {}

    points = tracer.output.points.to_array()

    lines = tracer.output.lines.to_array()

    for i in range(tracer.output.cell_data.number_of_arrays):
        arr = tracer.output.cell_data.get_array(i).to_array()
        name = tracer.output.cell_data.get_array_name(i)
        data[name] = arr
    start = 0
    line_segments = []
    point_idx = []
    for i in range(tracer.output.lines.number_of_cells):
        """loop over al lines"""
        n = lines[start]
        idx = lines[(start+1):(start+n+1)]
        line = points[idx]
        line_segments.append(line)
        point_idx.append(idx)
        start += (n + 1)
    data['line'] = line_segments
    data['points'] = point_idx
    return pd.DataFrame(data)

def make_particles(polydata, n=100, mode='random'):
    points = polydata.points.to_array()
    xmax, ymax, zmax = points.max(axis=0)
    xmin, ymin, zmin = points.min(axis=0)
    # take the srt of n
    n = np.round(np.sqrt(n)).astype('int')
    if mode != 'random':
        seed_X, seed_Y = np.meshgrid(
            np.linspace(xmin, xmax, num=n),
            np.linspace(ymin, ymax, num=n)
        )
    else:
        seed_X = np.random.random((n, n)) * (xmax - xmin) + xmin
        seed_Y = np.random.random((n, n)) * (ymax - ymin) + ymin
    seed_Z = np.zeros_like(seed_X)
    seed_points = np.c_[seed_X.ravel(), seed_Y.ravel(), seed_Z.ravel()]
    seed = tvtk.PolyData()
    seed.points = seed_points
    return seed

def export_lines(lines, filename):

    features = []
    for i, line in lines.iterrows():
        # convert to wgs84 or something ...
        lon, lat = line['line'][:,0], line['line'][:,1]
        linestring = geojson.LineString(coordinates=np.c_[lon, lat].tolist())
        properties = dict(line[['SeedIds', 'ReasonForTermination']])
        feature = geojson.Feature(id=line['SeedIds'], geometry=linestring, properties=properties)
        features.append(feature)
    fc = geojson.FeatureCollection(features=features)
    with open(filename, 'w') as f:
        geojson.dump(fc, f)
