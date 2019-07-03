from sympy import Array, diag, eye, simplify, sin, symbols
from sympy.tensor.tensor import TensExpr, TensMul

from einsteinpy.symbolic.metric import *
from einsteinpy.symbolic.tensor import *


def _generate_variables_simple():
    coords = symbols("x y z", real=True)
    # do not use any tensor related attributes of metric.
    # there is a dependency cycle between Tensor and Metric.
    metric = Metric("metric", coords, eye(3))
    return (coords, metric)


def _generate_schwarzschild():
    coords = symbols("t r theta phi", real=True)
    t, r, th, ph = coords
    schw = diag(1 - 1 / r, -1 / (1 - 1 / r), -r ** 2, -r ** 2 * sin(th) ** 2)
    g = Metric("g", coords, schw)
    mu, nu = indices("mu nu", g)
    return (coords, t, r, th, ph, schw, g, mu, nu)


def _generate_minkowski():
    coords = symbols("t x y z", real=True)
    t, x, y, z = coords
    mink = diag(1, -1, -1, -1)
    eta = Metric("eta", coords, mink)
    mu, nu = indices("mu nu", eta)
    return (coords, t, x, y, z, mink, eta, mu, nu)


def test_Index():
    (coords, metric) = _generate_variables_simple()
    mu = Index("mu", metric)
    assert isinstance(mu, TensorIndex)
    assert mu.is_up
    assert not (-mu).is_up


def test_indices():
    (coords, metric) = _generate_variables_simple()
    idxs = indices("i0:4", metric)
    assert len(idxs) == 4
    assert all([isinstance(i, Index) for i in idxs])


def test_Tensor():
    from sympy.tensor.tensor import tensorsymmetry

    (coords, metric) = _generate_variables_simple()
    T = Tensor("T", coords, metric)
    assert isinstance(T, TensorHead)
    assert T.as_array() == Array(coords)
    assert T.covar == (1,)
    assert T.symmetry == tensorsymmetry([1])


def test_Tensor_comm():
    (coords, metric) = _generate_variables_simple()
    T = Tensor("T", coords, metric)
    assert T.commutes_with(T) == 0
    assert T.commutes_with(metric) == 0
    assert T.commutes_with(metric.partial) is None


def test_Tensor_simplify():
    (coords, metric) = _generate_variables_simple()
    x, y, z = coords
    expr = -(1 - x) ** 2 / (x - 1)
    expr_simplified = 1 - x
    T = Tensor("T", 3 * [expr], metric)
    assert str(T.simplify()[0]) == str(expr_simplified)


def test_Tensor_covariance_transform():
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    E, p1, p2, p3 = symbols("E p_1:4", positive=True)
    p = Tensor("p", [E, p1, p2, p3], g)
    res = p.covariance_transform(-mu)
    expect = Array(
        [E * (1 - 1 / r), -p1 / (1 - 1 / r), -p2 * r ** 2, -p3 * r ** 2 * sin(th) ** 2]
    )
    assert res[0].equals(expect[0])
    assert res[1].equals(expect[1])
    assert res[2].equals(expect[2])
    assert res[3].equals(expect[3])


def test_AbstractTensor():
    (coords, metric) = _generate_variables_simple()
    T = Tensor("T", coords, metric)
    assert isinstance(T, AbstractTensor)
    assert isinstance(T.as_array(), Array)
    assert T.is_Tensor
    assert not T.is_Metric
    assert not T.is_Spacetime


def test_IndexedTensor():
    (coords, metric) = _generate_variables_simple()
    T = Tensor("T", coords, metric)
    mu = Index("mu", metric)
    assert isinstance(T(mu), IndexedTensor)


def test_ReplacementManager():
    (coords, metric) = _generate_variables_simple()
    T = Tensor("T", coords, metric)
    mu = Index("mu", metric)
    assert ReplacementManager.has(T)
    assert ReplacementManager.get_value(T(mu)) == T.as_array()


def test_expand_tensor():
    (coords, t, r, th, ph, mink, eta, mu, nu) = _generate_minkowski()
    E1, E2, E3, B1, B2, B3 = symbols("E_1:4 B_1:4", real=True)
    matrix = [[0, -E1, -E2, -E3], [E1, 0, -B3, B2], [E2, B3, 0, -B1], [E3, -B2, B1, 0]]
    F = Tensor("F", matrix, eta, symmetry=[[2]])
    assert expand_tensor(eta(mu, nu) * eta(-mu, -nu)) == 4
    assert expand_tensor(F(mu, -mu)) == 0
    assert (
        expand_tensor(F(mu, nu) * F(-mu, -nu))
        == 2 * B1 ** 2
        + 2 * B2 ** 2
        + 2 * B3 ** 2
        - 2 * E1 ** 2
        - 2 * E2 ** 2
        - 2 * E3 ** 2
    )
    (coords, t, r, th, ph, schw, g, mu, nu) = _generate_schwarzschild()
    x = Tensor("x", [t, r, th, ph], g)
    res = expand_tensor(g(mu, nu))
    assert schw.inv().equals(res)
    res1 = expand_tensor(g(mu, nu) * g(-mu, -nu))
    res2 = expand_tensor(g(-mu, -nu) * g(mu, nu))
    assert res1 == res2
    assert simplify(res1) == 4
    res1 = expand_tensor(g(-mu, -nu) * x(nu))
    res2 = expand_tensor(x(-mu))
    res3 = x.covariance_transform(-mu)
    assert res1 == res2
    assert res2 == res3
    res1 = expand_tensor(g(mu, nu) * x(-nu))
    res2 = expand_tensor(x(mu))
    assert simplify(res1) == res2
    res = expand_tensor(x(mu) * x(-mu))
    assert (
        res
        == t ** 2 * (1 - 1 / r)
        - r ** 2 / (1 - 1 / r)
        - th ** 2 * r ** 2
        - ph ** 2 * r ** 2 * sin(th) ** 2
    )
