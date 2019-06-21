from sympy import Array, cos, sin, symbols

from einsteinpy.symbolic.tensor import Tensor


def test_Tensor():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    test_arr = Array(test_list)
    obj1 = Tensor(test_arr)
    obj2 = Tensor(test_list)
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
