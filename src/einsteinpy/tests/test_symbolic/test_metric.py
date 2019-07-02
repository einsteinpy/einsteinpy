from sympy import Matrix, diag, symbols, sin, simplify
from sympy.tensor.tensor import TensMul
from einsteinpy.symbolic.tensor import *
from einsteinpy.symbolic.metric import *


def test_Metric():
    g = SpacetimeMetric("g", _coords, _schw)
    mu, nu = indices("mu nu", g)
    assert isinstance(g(mu, nu), IndexedTensor)
    res1, res2 = map(simplify, (g.covariance_transform(-mu, -nu), g.as_inverse()))
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
