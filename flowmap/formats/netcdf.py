import logging
import hashlib
import uuid
import json
import pathlib

import osgeo.osr
import mako.template

logger = logging.getLogger(__name__)


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


def file2uuid(fname):
    """read a file and return the unique uuid based on the file content"""
    hash_md5 = hashlib.md5()
    logger.info("computing uuid based on md5sum")
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    digest = hash_md5.digest()
    uuid_ = uuid.UUID(bytes=digest)
    logger.info("uuid: %s", uuid_)

    return uuid_


class NetCDF(object):
    def __init__(self, path, src_epsg=4326, dst_epsg=28992, vmin=-0.5, vmax=0.5, framescale=3.0):
        self.path = path
        # source and destination epsg code
        self.src_epsg = src_epsg
        self.dst_epsg = dst_epsg
        self.vmin = vmin
        self.vmax = vmax
        self.framescale = framescale
        logger.debug("Object constructed with %s", vars(self))

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

    def meta(self):
        metadata = {
            'metadata': {}
        }

        default_path = pathlib.Path('defaults.json')
        if default_path.exists():
            meta = json.load(default_path.open())
            metadata.update(meta)
        id_ = file2uuid(self.path)
        metadata['id'] = id_
        return metadata
