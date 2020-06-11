import numpy as np
import sympy

from einsteinpy.symbolic import ChristoffelSymbols, MetricTensor, RiemannCurvatureTensor


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


def test_circular_index_config_same():
    sch = schwarzschild_metric()
    rm0 = RiemannCurvatureTensor.from_metric(sch)
    rm1 = rm0.change_config("lulu")
    rm2 = rm1.change_config("uuuu")
    rm3 = rm2.change_config("ulll")
    dims = sch.dims
    for i in range(dims ** 4):
        # t,s,r,n each goes from 0 to (dims-1)
        # hack for codeclimate. Could be done with 4 nested for loops
        n = i % dims
        r = (int(i / dims)) % (dims)
        s = (int(i / (dims ** 2))) % (dims)
        t = (int(i / (dims ** 3))) % (dims)
        assert sympy.simplify(rm3[n, r, s, t] - rm0[n, r, s, t]) == 0
    assert rm0.parent_metric == rm0._parent_metric


def test_no_exception_raised_passing_christoffel_with_unconventional_indices():
    sch = schwarzschild_metric()
    ch = ChristoffelSymbols.from_metric(sch)
    ch2 = ch.change_config("uuu")
    try:
        rm = RiemannCurvatureTensor.from_christoffels(ch2)
        assert True
    except:
        assert False
