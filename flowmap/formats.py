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

    def variables(self):
        with netCDF4.Datset(self.path) as ds:
            variables = ds.variables.keys()
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
