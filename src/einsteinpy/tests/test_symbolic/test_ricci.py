import numpy as np
import sympy

from einsteinpy.symbolic import MetricTensor, RicciTensor


def spherical_metric():
    symbolstr = "r theta phi"
    syms = sympy.symbols(symbolstr)
    list2d = (np.zeros(shape=(3, 3), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sympy.sin(syms[1]) ** 2)
    metric = MetricTensor(list2d, syms)
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
