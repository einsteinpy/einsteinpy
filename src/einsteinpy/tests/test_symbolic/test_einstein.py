import numpy as np
import sympy
from sympy import cos, sin, sinh
from sympy.diffgeom import CovarDerivativeOp

from einsteinpy.symbolic import EinsteinTensor, MetricTensor


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


def test_einstein_ValueError_wrong_config_length():
    testarr = np.ones((4, 4), dtype=int).tolist()
    syms = sympy.symbols("x y")
    try:
        obj = EinsteinTensor(testarr, syms, config="lll", parent_metric=None)
        assert False
    except ValueError:
        assert True
