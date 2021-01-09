import sys
import types
import warnings
from importlib import reload
from unittest import mock

import einsteinpy.ijit


@mock.patch.dict(sys.modules, {"numba": None})
def test_warning_and_returntype():
    with warnings.catch_warnings(record=True) as w:
        custom_jit = reload(einsteinpy.ijit)

        @custom_jit.jit
        def _simple_func(a, b):
            return a + b

        assert len(w) >= 1
        assert isinstance(_simple_func, types.FunctionType)

@mock.patch.dict(sys.modules, {"numba": None})
def test_decorator():
    custom_jit = reload(einsteinpy.ijit)

    def _simple_func(a, b):
        return a + b

    res = custom_jit.jit()

    assert res.__name__ == "_jit"
    assert res(_simple_func).__name__ == "_simple_func"
