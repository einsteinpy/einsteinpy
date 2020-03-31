import numpy as np
import pytest
from sympy import Array, ImmutableDenseNDimArray, cos, sin, symbols
from sympy.abc import p, q, x, y

from einsteinpy.symbolic import TransformationMatrix, simplify_sympy_array
from einsteinpy.symbolic.helpers import _change_name

# tests for TransformationMatrix


@pytest.mark.parametrize(
    "iterable, old_coords, new_coords",
    [
        [[[x, x * y], [1, 2 * y]], (p, q), (x, y)],
        pytest.param([1, 2 * y], (p, q), (x, y), marks=pytest.mark.xfail),
    ],
)
def test_TransformationMatrix_is_initializing(iterable, old_coords, new_coords):
    obj = TransformationMatrix(iterable, old_coords, new_coords)
    assert isinstance(obj, ImmutableDenseNDimArray)
    assert obj.rank() == 2 and isinstance(obj._array, list)


@pytest.fixture
def cart2sph_matrix():
    old_coords, new_coords = [x, y], [p, q]
    new2old = [p * cos(q), p * sin(q)]
    obj = TransformationMatrix.from_new2old(old_coords, new_coords, new2old)
    return obj


def test_TransformationMatrix_classmethod_new2old(cart2sph_matrix):
    obj = cart2sph_matrix

    cmp_array = Array([[1 / cos(q), 1 / sin(q)], [-1 / (p * sin(q)), 1 / (p * cos(q))]])
    zero_array = Array(np.zeros(shape=cmp_array.shape, dtype=int))

    assert simplify_sympy_array(cmp_array - obj) == zero_array


def test_TransformationMatrix_inv_function(cart2sph_matrix):
    assert cart2sph_matrix._inv is None

    obj = cart2sph_matrix.inv()
    assert cart2sph_matrix._inv is not None

    obj = cart2sph_matrix.inv()

    cmp_array = Array([[cos(q), sin(q)], [-1 * p * sin(q), p * cos(q)]])
    zero_array = Array(np.zeros(shape=cmp_array.shape, dtype=int))

    assert simplify_sympy_array(cmp_array - obj) == zero_array


def test_TransformationMatrix_variables_are_not_None(cart2sph_matrix):
    # this test was required to mitogate confusion arised
    # due to usage of both __new__ and __init__
    obj = cart2sph_matrix
    assert (obj.old_coords is not None) and (obj.new_coords is not None)
    assert obj.new2old is not None


@pytest.mark.parametrize(
    "curr_name, context, expected",
    [
        ["tens", "__lt", "tens__lt"],
        ["tens__ulu", "__ull", "tens__ulu__ull"],
        ["tens__lt", "__lt", "tens__lt__lt"],
        [None, "__lt", None],
    ],
)
def test_change_name(curr_name, context, expected):
    altered_name = _change_name(curr_name, context)
    assert altered_name == expected
