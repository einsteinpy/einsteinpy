from sympy import Expr, diag, sin, symbols, zeros

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


def test_Metric():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    assert isinstance(g, AbstractTensor)
    assert isinstance(g, TensorIndexType)
    assert isinstance(g(mu, nu), IndexedTensor)
    assert isinstance(g.metric, Tensor)
    assert isinstance(g.partial, PartialDerivative)
    assert g.metric.covar == (-1, -1)


def test_SpacetimeMetric():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    array1 = g.as_array()
    assert g.signature == (1, -1, -1, -1)
    rev_sig = g.reverse_signature()
    array2 = g.as_array()
    assert rev_sig == (-1, 1, 1, 1)
    assert array1 == -1 * array2


def test_Metric_christoffel():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    gamma = g.christoffel.simplify()
    assert (
        str(gamma)
        == "[[[0, 0.5/(r*(r - 1)), 0, 0], [0.5/(r*(r - 1)), 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[0.5*(r - 1)/r**3, 0, 0, 0], [0, -0.5/(r*(r - 1)), 0, 0], [0, 0, 1.0 - 1.0*r, 0], [0, 0, 0, (1.0 - 1.0*r)*sin(theta)**2]], [[0, 0, 0, 0], [0, 0, 1.0/r, 0], [0, 1.0/r, 0, 0], [0, 0, 0, -0.5*sin(2*theta)]], [[0, 0, 0, 0], [0, 0, 0, 1.0/r], [0, 0, 0, 1.0/tan(theta)], [0, 1.0/r, 1.0/tan(theta), 0]]]"
    )


def test_Metric_riemann():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    rh, si = indices("rho sigma", g)
    R = g.riemann

    def is_zero(arr):
        for comp in arr:
            yield comp.equals(0)

    expr = R(-rh, -si, -mu, -nu) + R(-si, -rh, -mu, -nu)
    assert all(is_zero(expand_tensor(expr)))
    expr = R(-rh, -si, -mu, -nu) + R(-rh, -si, -nu, -mu)
    assert all(is_zero(expand_tensor(expr)))
    expr = R(-rh, -si, -mu, -nu) - R(-mu, -nu, -rh, -si)
    assert all(is_zero(expand_tensor(expr)))
    expr = R(-rh, -si, -mu, -nu) + R(-rh, -mu, -nu, -si) + R(-rh, -nu, -si, -mu)
    assert all(is_zero(expand_tensor(expr)))


def test_Metric_ricci_tensor():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    R = g.ricci_tensor
    assert zeros(4).equals(R.as_array())


def test_Metric_ricci_scalar():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    R = g.ricci_scalar
    assert isinstance(R, Expr)
    assert R.equals(0)


def test_Metric_einstein():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    G = g.einstein
    assert zeros(4).equals(G.as_array())


def test_Metric_weyl():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    rh, si = indices("rho sigma", g)
    C = g.weyl

    def is_zero(arr):
        for comp in arr:
            yield comp.equals(0)

    expr = C(-rh, -si, -mu, -nu) + C(-si, -rh, -mu, -nu)
    assert all(is_zero(expand_tensor(expr)))
    expr = C(-rh, -si, -mu, -nu) + C(-rh, -si, -nu, -mu)
    assert all(is_zero(expand_tensor(expr)))
    expr = C(-rh, -si, -mu, -nu) - C(-mu, -nu, -rh, -si)
    assert all(is_zero(expand_tensor(expr)))
    expr = C(-rh, -si, -mu, -nu) + C(-rh, -mu, -nu, -si) + C(-rh, -nu, -si, -mu)
    assert all(is_zero(expand_tensor(expr)))
