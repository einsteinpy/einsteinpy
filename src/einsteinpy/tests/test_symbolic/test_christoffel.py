import numpy as np
import sympy

from einsteinpy.symbolic import ChristoffelSymbols, MetricTensor


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


def test_change_config():
    x, y, z = sympy.symbols("x y z")

    list3d = np.zeros((3, 3, 3), dtype=int).tolist()
    for i in range(3):
        list3d[i][i][i] = (x ** i) * (y * (2 - i)) + i * z
    list3d[1][2][0] = list3d[1][0][2] = x * y * z
    list3d[2][1][0] = list3d[2][0][1] = 4 * y

    metriclist = np.identity(3).tolist()
    metric = MetricTensor(metriclist, (x, y, z), "uu")
    ch = ChristoffelSymbols(list3d, (x, y, z), "ull", parent_metric=metric)
    chr_new = ch.change_config("llu")

    for t in range(3):
        i, j, k = t % 3, (int(t / 3)) % 3, (int(t / (3 ** 2))) % 3
        assert sympy.simplify(ch[i, j, k] - chr_new[i, j, k]) == 0


def test_wrong_number_of_indices_ValueError():
    x, y, z = sympy.symbols("x y z")

    list3d = np.zeros((3, 3, 3), dtype=int).tolist()
    for i in range(3):
        list3d[i][i][i] = (x ** i) * (y * (2 - i)) + i * z
    list3d[1][2][0] = list3d[1][0][2] = x * y * z
    list3d[2][1][0] = list3d[2][0][1] = 4 * y

    try:
        ch = ChristoffelSymbols(list3d, (x, y, z), "ulll")
        assert False
    except ValueError:
        assert True


def test_properties():
    sch_inv = schwarzschild_metric()
    ch = ChristoffelSymbols.from_metric(sch_inv)
    assert ch.parent_metric == ch._parent_metric
    # test change_config, should raise ValueError
    ch._parent_metric = None
    try:
        ch_new = ch.change_config("lll")
        assert False
    except Exception:
        return True
