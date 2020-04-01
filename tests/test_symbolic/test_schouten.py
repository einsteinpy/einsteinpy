import numpy as np
import sympy
from sympy import cos, sin, sinh

from einsteinpy.symbolic import (
    MetricTensor,
    RicciScalar,
    RicciTensor,
    RiemannCurvatureTensor,
    SchoutenTensor,
    simplify_sympy_array,
)


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


def schwarzschild_metric():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    G, M, c, a = sympy.symbols("G M c a")
    # using metric values of schwarschild space-time
    # a is schwarzschild radius
    list2d = np.zeros((4, 4), dtype=int).tolist()
    list2d[0][0] = 1 - (a / syms[1])
    list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
    list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
    list2d[3][3] = -1 * (syms[1] ** 2) * (sympy.sin(syms[2]) ** 2) / (c ** 2)
    sch = MetricTensor(list2d, syms)
    # print (sch.tensor())
    return sch


def test_SchoutenTensor_parent_metric_property_and_ValueError():
    mt = anti_de_sitter_metric()
    St = SchoutenTensor.from_metric(mt)
    assert St.parent_metric == St._parent_metric and St.parent_metric == mt
    # ValueError test part
    try:
        St2 = SchoutenTensor(St.tensor(), St.syms, config="uuu")
        assert False
    except ValueError:
        assert True


def test_Schouten_dim2():
    list_metric = np.zeros((2, 2), dtype=int).tolist()
    symbolstr = "t r"
    syms = sympy.symbols(symbolstr)
    metric = MetricTensor(list_metric, syms)
    try:
        obj = SchoutenTensor.from_metric(metric)
        assert False
    except ValueError:
        assert True


def test_SchoutenTensor_with_arbritary_config():
    sch = schwarzschild_metric()
    ch = SchoutenTensor.from_metric(sch)
    assert ch.parent_metric == ch._parent_metric
    # test change_config, should raise ValueError
    ch._parent_metric = None
    try:
        ch_new = ch.change_config("lll")
        assert False
    except Exception:
        return True


def test_SchoutenTensor_config_change():
    metric = schwarzschild_metric()
    t1 = SchoutenTensor.from_metric(metric)
    t2 = t1.change_config("ul")
    t3 = t2.change_config("ll")
    assert t1.config == "ll" and t3.config == "ll"
    compare = sympy.Array(np.zeros((4, 4), dtype=int))
    testarr = simplify_sympy_array(t1.tensor() - t3.tensor())
    assert testarr == compare
