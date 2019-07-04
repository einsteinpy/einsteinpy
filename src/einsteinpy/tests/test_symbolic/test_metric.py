import numpy as np
import sympy

from einsteinpy.symbolic import MetricTensor


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


def test_MetricTensor():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    obj = schwarzschild_metric()
    assert obj.dims == 4
    assert obj.symbols() == syms


def test_TypeError():
    list2d = np.zeros((4, 4), dtype=int).tolist()
    syms = 100
    try:
        obj = MetricTensor(list2d, syms)
        assert False
    except TypeError:
        assert True


def test_properties():
    sch = schwarzschild_metric()
    assert sch.config == "ll"
    assert sch.order == 2


def test_inv():
    sch = schwarzschild_metric()
    sch_inv = sch.inv()
    assert sch == sch_inv.inv()
    assert sch_inv == sch.inv()
    for i in range(4):
        # this is because schwarzschild metric is diagonal
        assert sympy.simplify(sch[i, i] * sch_inv[i, i]) == 1


def test_change_config_ValueError():
    sch = schwarzschild_metric()
    try:
        sch_new = sch.change_config("ul")
        assert False
    except ValueError:
        assert True


def test_check_ValueError_on_incorrect_config_length():
    sch = schwarzschild_metric()
    try:
        sch_new = MetricTensor(sch.tensor(), sch.syms, config="lll")
        assert False
    except ValueError:
        assert True


def test_use_original_config_in_change_config():
    sch = schwarzschild_metric()
    sch_new = sch.change_config("ll")
    assert sch == sch_new
