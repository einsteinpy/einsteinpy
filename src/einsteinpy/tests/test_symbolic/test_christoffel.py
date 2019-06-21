import numpy as np
import sympy

from einsteinpy.symbolic import ChristoffelSymbols, MetricTensor


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
    # print (sch.tensor())
    return sch


def test_ChristoffelSymbols():
    sch = schwarzschild_metric()
    chl = ChristoffelSymbols.from_metric(sch)
    mat = chl.tensor()
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    G, M, c, a = sympy.symbols("G M c a")
    assert (
        sympy.simplify(mat[2, 3, 3] - (-1 * sympy.cos(syms[2]) * sympy.sin(syms[2])))
        == 0
    )
    assert sympy.simplify(mat[3, 3, 1] - syms[1] / (syms[1] ** 2)) == 0
    assert (
        sympy.simplify(
            (mat[1, 1, 1].subs({a: (2 * G * M / (c ** 2))}))
            - (G * M / (2 * G * M * syms[1] - c ** 2 * syms[1] ** 2))
        )
        == 0
    )
    assert chl.symbols() == syms


def test_TypeError():
    testarr = np.ones((4, 4, 4), dtype=int).tolist()
    syms = 0
    try:
        obj = ChristoffelSymbols(testarr, syms)
        assert False
    except TypeError:
        assert True
