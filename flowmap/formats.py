# supported file formats
import netCDF4
import logging

logger = logging.getLogger(__name__)

# what to export on import from *, also used to get a list of available formats
__all__ = ["Matroos"]


class Matroos(object):
    """FEWS Matroos format"""
    def __init__(self, path):
        self.path = path

    def globals(self):
        """generate global variables"""
        with netCDF4.Datset(self.path) as ds:
            times = netCDF4.num2date(ds.variables['time'][:],
                                     ds.variables['time'].units)
            analysis_time = netCDF4.num2date(ds.variables['analysis_time'][:],
                                             ds.variables['analysis_time'].units)

            lat = ds.variables['lat'][:]
            lon = ds.variables['lon'][:]

            # initial values (used to determine shapes and stuff, maybe remove if not used)
            sep = ds.variables['sep'][0]
            u1 = ds.variables['velu'][0]
            v1 = ds.variables['velv'][0]

        variables = dict(
            times=times,
            analysis_time=analysis_time,
            u1=u1,
            v1=v1,
            sep=sep,
            lat=lat,
            lon=lon
        )
        return variables

    def validate(self):
        """validate a file"""
        valid = True
        with netCDF4.Datset(self.path) as ds:
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
