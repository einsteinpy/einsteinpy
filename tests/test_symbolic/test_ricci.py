import numpy as np
import sympy
from sympy import cos, sin, sinh

from einsteinpy.symbolic import (
    MetricTensor,
    RicciScalar,
    RicciTensor,
    RiemannCurvatureTensor,
)


def spherical_metric():
    symbolstr = "r theta phi"
    syms = sympy.symbols(symbolstr)
    list2d = (np.zeros(shape=(3, 3), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sympy.sin(syms[1]) ** 2)
    metric = MetricTensor(list2d, syms)
    return metric


def anti_de_sitter_metric():
    coords = sympy.symbols("t chi theta phi")
    t, ch, th, ph = coords
    m = sympy.diag(
        -1,
        cos(t) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2 * sin(th) ** 2,
    ).tolist()
    metric = MetricTensor(m, coords)
    return metric


def test_Ricci_is_zero_for_spherical_space():
    metric = spherical_metric()
    dims = metric.dims
    T = RicciTensor.from_metric(metric)
    arr = T.tensor()
    for i in range(dims ** 2):
        n = i % dims
        r = (int(i / dims)) % (dims)
        assert sympy.simplify(arr[n, r]) == 0


def test_TypeError():
    testarr = np.ones((4, 4), dtype=int).tolist()
    syms = 0
    try:
        obj = RicciTensor(testarr, syms)
        assert False
    except TypeError:
        assert True


def test_RicciTensor_parent_metric_property_and_ValueError():
    mt = anti_de_sitter_metric()
    Rt = RicciTensor.from_metric(mt)
    assert Rt.parent_metric == Rt._parent_metric and Rt.parent_metric == mt
    # ValueError test part
    try:
        Rt2 = RicciTensor(Rt.tensor(), Rt.syms, config="uuu")
        assert False
    except ValueError:
        assert True


def test_RicciTensor_from_riemann_with_arbritary_config():
    mt = anti_de_sitter_metric()
    rm = RiemannCurvatureTensor.from_metric(mt)
    rm2 = rm.change_config("lulu")
    try:
        Rt = RicciTensor.from_riemann(rm2)
        assert True
    except Exception:
        assert False


# Tests for Ricci Scalar


def test_RicciScalar_is_constant_for_ADS_spacetime():
    # R is -12 for anti-de-Sitter spacetime
    mt = anti_de_sitter_metric()
    R = RicciScalar.from_metric(mt)
    assert sympy.simplify(R.expr).is_constant()


def test_RicciScalar_raise_TypeError():
    x, y = sympy.symbols("x y")
    try:
        R = RicciScalar(x ** 2 * y * 3, 0)
        assert False
    except TypeError:
        assert True


def test_RicciScalar_inbuilt_simplify_works():
    x, y = sympy.symbols("x y")
    R = RicciScalar((-y) + 2 * sin(x) * cos(x) + y - sin(2 * x), (x, y))
    val = R.simplify(set_self=True)
    assert R.expr == 0
    assert sum(val) == 0
