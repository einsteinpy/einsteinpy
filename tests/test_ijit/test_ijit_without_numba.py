import sys
import types
import warnings
from importlib import reload
from unittest import mock

import einsteinpy.ijit


@mock.patch.dict(sys.modules, {"numba": None})
def test():
    with warnings.catch_warnings(record=True) as w:
        custom_jit = reload(einsteinpy.ijit)

        @custom_jit.jit
        def _simple_func(a, b):
            return a + b

        assert len(w) >= 1
        assert isinstance(_simple_func, types.FunctionType)


custom_jit = reload(einsteinpy.ijit)
