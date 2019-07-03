from sympy import Matrix, diag, simplify, sin, symbols
from sympy.tensor.tensor import TensMul

from einsteinpy.symbolic.metric import *
from einsteinpy.symbolic.tensor import *

_t, _r, _th, _ph = symbols("t r theta phi", real=True)
_coords = [_t, _r, _th, _ph]
_schw = diag(1 - 1 / _r, -1 / (1 - 1 / _r), -_r ** 2, -_r ** 2 * sin(_th) ** 2)
_E, _p1, _p2, _p3 = symbols("E p_1:4", positive=True)
_momentum = [_E, _p1, _p2, _p3]


def test_Metric():
    g = SpacetimeMetric("g", _coords, _schw)
    mu, nu = indices("mu nu", g)
    assert isinstance(g(mu, nu), IndexedTensor)
    res1, res2 = map(simplify, (g.covariance_transform(mu, nu), g.as_inverse()))
    assert res1 == res2


def test_SpacetimeMetric():
    g = SpacetimeMetric("g", _coords, _schw)
    array1 = g.as_array()
    assert g.signature == (1, -1, -1, -1)
    rev_sig = g.reverse_signature()
    array2 = g.as_array()
    assert rev_sig == (-1, 1, 1, 1)
    assert array1 == -1 * array2


def test_TypeError():
    pass
