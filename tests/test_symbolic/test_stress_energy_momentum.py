import numpy as np
import sympy
from sympy import cos, sin, sinh

from einsteinpy.symbolic import (
    MetricTensor,
    StressEnergyMomentumTensor,
    simplify_sympy_array,
)


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
    return sch


def test_StressEnergyMomentumTensor_ValueError_wrong_config_length():
    testarr = np.ones((4, 4), dtype=int).tolist()
    syms = sympy.symbols("t x y z")
    try:
        obj = StressEnergyMomentumTensor(
            testarr, syms, config="lll", parent_metric=None
        )
        assert False
    except ValueError:
        assert True


def test_StressEnergyMomentumTensor_change_config_no_parent_metric():
    testarr = np.ones((4, 4), dtype=int).tolist()
    syms = sympy.symbols("t x y z")
    obj = StressEnergyMomentumTensor(testarr, syms, config="ll")
    try:
        obj2 = obj.change_config("lu")
        assert False
    except Exception:
        assert True


def test_StressEnergyMomentumTensor_cyclic_config_change():
    metric = schwarzschild_metric()
    t1 = StressEnergyMomentumTensor.from_metric(metric)
    t2 = t1.change_config("ul")
    t3 = t2.change_config("ll")
    assert t1.config == "ll" and t3.config == "ll"
    compare = sympy.Array(np.zeros((4, 4), dtype=int))
    testarr = simplify_sympy_array(t1.tensor() - t3.tensor())
    assert testarr == compare
