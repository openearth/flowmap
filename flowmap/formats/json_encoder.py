# -*- coding: utf-8 -*-
import json

registry = []

try:
    import numpy as np

    def numpy_default(obj):
        """Support for data types that JSON default encoder
            does not do.

        This includes:

            * Numpy array or number
            * Complex number
            * Set
            * Bytes (Python 3)

        Examples
        --------
        >>> import json
        >>> import numpy as np
        >>> from astropy.utils.misc import JsonCustomEncoder
        >>> json.dumps(np.arange(3), cls=JsonCustomEncoder)
        '[0, 1, 2]'

        """
        if isinstance(obj, (np.ndarray, np.number)):
            return obj.tolist()
        elif isinstance(obj, (complex, np.complex)):
            return [obj.real, obj.imag]
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, bytes):  # pragma: py3
            return obj.decode()
        else:
            raise TypeError("not supported")
    registry.append(numpy_default)
except ImportError:
    pass

try:
    import pandas as pd

    def pandas_default(obj):
        try:
            return json.loads(obj.to_json(orient='records'))
        except AttributeError:
            raise TypeError("not a to_json function")
    registry.append(pandas_default)
except ImportError:
    pass


class RegistryEncoder(json.JSONEncoder):
    def default(self, obj):
        # try with all default functions in registry
        for default in registry:
            try:
                result = default(obj)
                return result
            except TypeError as e:
                pass
        return super(json.JSONEncoder, self).default(obj)
