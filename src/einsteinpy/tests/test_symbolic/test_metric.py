import sympy

from einsteinpy.symbolic import MetricTensor


def schwarzschild_metric():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    G, M, c, a = sympy.symbols("G M c a")
    # using metric values of schwarschild space-time
    # a is schwarzschild radius
    list2d = [[0 for j in range(4)] for i in range(4)]
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
    list2d = [[0 for i in range(4)] for i in range(4)]
    syms = 100
    try:
        obj = MetricTensor(list2d, syms)
        assert False
    except TypeError:
        assert True
