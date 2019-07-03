from sympy import diag, sin, symbols

from einsteinpy.symbolic.metric import *
from einsteinpy.symbolic.partial import *
from einsteinpy.symbolic.tensor import *


def _generate_schwarzschild():
    coords = symbols("t r theta phi", real=True)
    t, r, th, ph = coords
    schw = diag(1 - 1 / r, -1 / (1 - 1 / r), -r ** 2, -r ** 2 * sin(th) ** 2)
    g = SpacetimeMetric("g", coords, schw, timelike=True)
    mu, nu = indices("mu nu", g)
    return (coords, t, r, th, ph, schw, g, mu, nu)


def test_DiffOperator():
    r = Symbol("r")
    dr = DiffOperator(r)
    expr = dr * (1 - 1 / r)
    assert expr == 1 / r ** 2
    expr = (1 - 1 / r) * dr
    assert isinstance(expr, DiffOperator)
    assert (expr * (r ** 3 / 3)).equals(r * (r - 1))


def test_PartialDerivative():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    d = g.partial
    assert isinstance(d, PartialDerivative)
    assert isinstance(d, Tensor)
    assert d.commutes_with(d) == 0
    assert d.covar == (-1,)
    assert all([isinstance(dx, DiffOperator) for dx in d.as_array()])
