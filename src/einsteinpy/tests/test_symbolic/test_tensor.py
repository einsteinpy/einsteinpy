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
    for i in range(len(test_list)):
        assert obj1[i]
    for i in range(len(test_arr)):
        assert obj2[i]
    assert obj1.tensor() == obj2.tensor()
    assert isinstance(obj1.tensor(), Array)
    assert obj1.simplify()[0, 1, 1] == 0


def test_TypeError():
    arr = 0
    try:
        obj = Tensor(arr)
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
