import numpy as np
from sympy import Array, cos, simplify, sin, symbols

from einsteinpy.symbolic.tensor import Tensor


def schwarzschild_tensor():
    symbolstr = "t r theta phi"
    syms = symbols(symbolstr)
    G, M, c, a = symbols("G M c a")
    # using metric values of schwarschild space-time
    # a is schwarzschild radius
    list2d = np.zeros((4, 4), dtype=int).tolist()
    list2d[0][0] = 1 - (a / syms[1])
    list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
    list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
    list2d[3][3] = -1 * (syms[1] ** 2) * (sin(syms[2]) ** 2) / (c ** 2)
    sch = Tensor(list2d)
    return sch


def test_Tensor():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    test_arr = Array(test_list)
    obj1 = Tensor(test_arr)
    obj2 = Tensor(test_list)
    assert obj1.tensor() == obj2.tensor()
    assert isinstance(obj1.tensor(), Array)
    assert obj1.simplify()[0, 1, 1] == 0


def test_Tensor_getitem():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    obj = Tensor(test_list)
    n = 2
    for i in range(n ** 3):
        p, q, r = i % n, int(i / n) % n, int(i / n ** 2) % n
        assert obj[p, q, r] - test_list[p][q][r] == 0


def test_TypeError1():
    arr = 0
    try:
        obj = Tensor(arr)
        assert False
    except TypeError:
        assert True


def test_TypeError2():
    scht = schwarzschild_tensor().tensor()
    # pass non str object
    try:
        obj = Tensor(scht, config=0)
        assert False
    except TypeError:
        assert True


def test_TypeError3():
    scht = schwarzschild_tensor().tensor()
    # pass string containing elements other than 'l' or 'u'
    try:
        obj = Tensor(scht, config="al")
        assert False
    except TypeError:
        assert True


def test_subs_single():
    # replacing only schwarzschild radius(a) with 0
    T = schwarzschild_tensor()
    a, c = symbols("a c")
    test_arr = T.subs(a, 0)
    assert simplify(test_arr.arr[0, 0] - 1) == 0
    assert simplify(test_arr.arr[1, 1] - (-1 / c ** 2)) == 0


def test_subs_multiple():
    # replacing a with 0, c with 1
    # this should give a metric for spherical coordinates
    T = schwarzschild_tensor()
    a, c = symbols("a c")
    test_arr = T.subs([(a, 0), (c, 1)])
    assert simplify(test_arr.arr[0, 0] - 1) == 0
    assert simplify(test_arr.arr[1, 1] - (-1)) == 0


def test_check_properties():
    T = schwarzschild_tensor()
    assert T.order == T._order
    assert T.config == T._config
