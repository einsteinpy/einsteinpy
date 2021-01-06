import pytest
import numba

from einsteinpy.ijit import ijit


@pytest.mark.parametrize("test_input, expected", [(None, "_jit"), (numba.njit, "njit")])
def test_correct_func_returned(test_input, expected):
    returned_func = ijit(test_input)
    assert returned_func.__name__ == expected
