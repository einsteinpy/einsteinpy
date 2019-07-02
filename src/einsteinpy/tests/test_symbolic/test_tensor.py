from sympy import Matrix, diag, symbols, sin, simplify
from sympy.tensor.tensor import TensMul
from einsteinpy.symbolic.tensor import *
from einsteinpy.symbolic.metric import *

_t, _r, _th, _ph = symbols('t r theta phi', real=True)
_coords = [_t, _r, _th, _ph]
_schw = diag(1-1/_r, -1/(1-1/_r), -_r**2, -_r**2*sin(_th)**2)
_E, _p1, _p2, _p3 = symbols('E p_1:4', positive=True)
_momentum = [_E, _p1, _p2, _p3]

def test_AbstractTensor():
    g = SpacetimeMetric('g', _coords, _schw)
    p = Tensor('p', _momentum, g)
    assert all([t.is_Tensor for t in (g, p)])
    assert g.is_Metric
    assert g.is_Spacetime

def test_IndexedTensor():
    g = SpacetimeMetric('g', _coords, _schw)
    p = Tensor('p', _momentum, g)
    mu, nu = indices('mu nu', g)
    expr1 = p(mu) * p(nu)
    assert isinstance(expr1, TensMul)
    assert expr1.rank == 2
    expr2 = p(mu) * p(-mu)
    assert expr2.rank == 0

def test_Tensor():
    g = SpacetimeMetric('g', _coords, _schw)
    p = Tensor('p', _momentum, g)
    mu, nu = indices('mu nu', g)
    assert p.rank == 1
    assert isinstance(p(mu), IndexedTensor)
    assert Matrix(p.covariance_transform(-mu)) \
        == _schw.inv() @ Matrix(p.as_array())

def test_Metric():
    g = SpacetimeMetric('g', _coords, _schw)
    mu, nu = indices('mu nu', g)
    assert isinstance(g(mu,nu), IndexedTensor)
    res1, res2 = map(simplify, (g.covariance_transform(-mu, -nu), g.as_inverse()))
    assert res1 == res2

def test_SpacetimeMetric():
    g = SpacetimeMetric('g', _coords, _schw)
    array1 = g.as_array()
    assert g.signature == (1, -1, -1, -1)
    rev_sig = g.reverse_signature()
    array2 = g.as_array()
    assert rev_sig == (-1, 1, 1, 1)
    assert array1 == -1*array2

def test_Index():
    g = SpacetimeMetric('g', _coords, _schw)
    mu, nu = indices('mu nu', g)
    assert mu.is_up
    assert mu != -mu
    assert not (-mu).is_up

def test_expand_tensor():
    g = SpacetimeMetric('g', _coords, _schw)
    p = Tensor('p', _momentum, g)
    mu, nu = indices('mu nu', g)
    expr1 = g(-mu,-nu) * p(nu)
    expr2 = p(-mu)
    res1, res2 = map(expand_tensor, (expr1, expr2))
    assert res1 == res2
    expr3 = g(-mu,-nu) * p(mu) * p(nu)
    expr4 = p(mu) * p(-mu)
    res3, res4 = map(expand_tensor, (expr3, expr4))
    assert res3 == res4
    assert res3 == \
        _E**2*(1-1/_r) - _p1**2/(1-1/_r) - _p2**2*_r**2 - _p3**2*_r**2*sin(_th)**2
    expr5 = g(-mu, -nu) * g(mu, nu)
    assert expand_tensor(expr5) == 4
    expr6 = g(mu, nu) * g(-mu, -nu)
    assert expand_tensor(expr6) == 4

def test_TypeError():
    pass
