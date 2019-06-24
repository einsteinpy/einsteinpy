import numpy as np
import sympy

from einsteinpy.symbolic import MetricTensor, RiemannCurvatureTensor


def spherical_metric():
    symbolstr = "r theta phi"
    syms = sympy.symbols(symbolstr)
    list2d = (np.zeros(shape=(3, 3), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sympy.sin(syms[1]) ** 2)
    metric = MetricTensor(list2d, syms)
    return metric


def test_Riemann_is_zero_for_spherical_space():
    metric = spherical_metric()
    dims = metric.dims
    T = RiemannCurvatureTensor.from_metric(metric)
    arr = T.tensor()
    for i in range(dims ** 4):
        n = i % dims
        r = (int(i / dims)) % (dims)
        s = (int(i / (dims ** 2))) % (dims)
        t = (int(i / (dims ** 3))) % (dims)
        assert arr[t, s, r, n] == 0


def test_TypeError():
    testarr = np.ones((4, 4, 4, 4), dtype=int).tolist()
    syms = 0
    try:
        obj = RiemannCurvatureTensor(testarr, syms)
        assert False
    except TypeError:
        assert True
