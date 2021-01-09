import sys
import types
import warnings
from importlib import reload
from unittest import mock

import einsteinpy.ijit


@mock.patch.dict(sys.modules, {"numba": None})
def test_decorator():
    with warnings.catch_warnings(record=True) as w:
        custom_jit = reload(einsteinpy.ijit)

        @custom_jit.jit
        def _simple_func(a, b):
            return a + b

        assert len(w) >= 1
        assert isinstance(_simple_func, types.FunctionType)

@mock.patch.dict(sys.modules, {"numba": None})
def test_function_given():
    with warnings.catch_warnings(record=True) as w:
        custom_jit = reload(einsteinpy.ijit)

        def _simple_func(a, b):
            return a + b

        return_func = custom_jit.jit(_simple_func)

        assert len(w) >= 1
        assert return_func.__name__ == "_simple_func"


@mock.patch.dict(sys.modules, {"numba": None})
def test_no_function_given():
    with warnings.catch_warnings(record=True) as w:
        custom_jit = reload(einsteinpy.ijit)
        returned_func = custom_jit.jit()

        assert len(w) >= 1
        assert returned_func.__name__ == "_jit"
