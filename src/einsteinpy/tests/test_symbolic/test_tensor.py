import numpy as np
import pytest
from numpy.testing import assert_allclose
from sympy import Array, Function, cos, simplify, sin, symbols
from sympy.abc import y, z

#Making an xfail marker to indicate that you expect a test to fail

xfail = pytest.mark.xfail

from einsteinpy.symbolic import (
    BaseRelativityTensor,
    MetricTensor,
    Tensor,
    simplify_sympy_array,
)

# Test for simplify_sympy_array()

def zero_expression():
    return 2 * sin(y * z) * cos(z * y) - sin(2 * z * y)

@pytest.mark.parametrize(
    "target",
    (Array([zero_expression(), 3]), Array(zero_expression()), zero_expression()),
)

# Marking unexpectedly failing test functions
@xfail    
def test_simplify_sympy_array_works_for_all(target):
    try:
        simplify_sympy_array(target)
        assert True
    except Exception:
        assert False

# Tests for Tensor and BaseRelativityTensor

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

def schwarzschild_metric():
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
    sch = MetricTensor(list2d, syms)
    return sch

def arbitrary_tensor1():
    symbolstr = "x0 x1 x2 x3"
    syms = symbols(symbolstr)
    a, c = symbols("a c")
    f1, f2, f3 = Function("f1")(a, syms[2]), Function("f2")(c), Function("f3")
    list2d = np.zeros((4, 4), dtype=int).tolist()
    list2d[0][0] = 1 - (a * f1 / syms[1])
    list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
    list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
    list2d[3][3] = -1 * (syms[1] ** 2) * (sin(syms[2]) ** 2) / (c ** 2)
    list2d[0][3] = list2d[3][0] = 5 * f2
    list2d[2][1] = list2d[1][2] = f3
    return BaseRelativityTensor(list2d, syms, config="ll"), [a, c], [f1, f2, f3]

@xfail
def test_Tensor():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    test_arr = Array(test_list)
    obj1 = Tensor(test_arr)
    obj2 = Tensor(test_list)
    assert obj1.tensor() == obj2.tensor()
    assert isinstance(obj1.tensor(), Array)

@xfail
def test_Tensor_simplify():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    obj = Tensor(test_list)
    
    # with set_self = False
    
    assert obj.simplify(set_self=False)[0, 1, 1] == 0
    assert not obj.tensor()[0, 1, 1] == 0
    
    # with set_self = True
    
    obj.simplify(set_self=True)
    assert obj.tensor()[0, 1, 1] == 0

@xfail
def test_Tensor_getitem():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    obj = Tensor(test_list)
    n = 2
    for i in range(n ** 3):
        p, q, r = i % n, int(i / n) % n, int(i / n ** 2) % n
        assert obj[p, q, r] - test_list[p][q][r] == 0


@xfail
def test_Tensor_str():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, x]], [[z, x], [y, z]]]
    obj1 = Tensor(test_list)

    assert "object at 0x" not in str(obj1)

@xfail
def test_Tensor_repr():
    x, y, z = symbols("x y z")
    test_list = [[[x, y], [y, sin(2 * z) - 2 * sin(z) * cos(z)]], [[z ** 2, x], [y, z]]]
    obj1 = Tensor(test_list)

    machine_representation = repr(obj1)
    assert not "object at 0x" in machine_representation

@xfail(raises=TypeError)
def test_TypeError2():
    scht = schwarzschild_tensor().tensor()
    
    # pass non str object
    
    try:
        obj = Tensor(scht, config=0)
        assert False
    except TypeError:
        assert True


@xfail(raises=TypeError)
def test_TypeError3():
    scht = schwarzschild_tensor().tensor()
    
    # pass string containing elements other than 'l' or 'u'
    
    try:
        obj = Tensor(scht, config="al")
        assert False
    except TypeError:
        assert True

@xfail
def test_subs_single():
    
    # replacing only schwarzschild radius(a) with 0
    
    T = schwarzschild_tensor()
    a, c = symbols("a c")
    test_arr = T.subs(a, 0)
    assert simplify(test_arr.arr[0, 0] - 1) == 0
    assert simplify(test_arr.arr[1, 1] - (-1 / c ** 2)) == 0

@xfail
def test_subs_multiple():
    
    # replacing a with 0, c with 1
    # this should give a metric for spherical coordinates
    
    T = schwarzschild_tensor()
    a, c = symbols("a c")
    test_arr = T.subs([(a, 0), (c, 1)])
    assert simplify(test_arr.arr[0, 0] - 1) == 0
    assert simplify(test_arr.arr[1, 1] - (-1)) == 0

@xfail
def test_check_properties():
    T = schwarzschild_tensor()
    assert T.order == T._order
    assert T.config == T._config

# Tests for BaseRelativityTensor

@xfail
def test_BaseRelativityTensor_automatic_calculation_of_free_variables():
    t1, variables, functions = arbitrary_tensor1()
    t2 = BaseRelativityTensor(
        t1.arr, t1.symbols(), config=t1.config, variables=variables, functions=functions
    )
    assert len(t1.variables) == len(t2.variables) and len(t1.variables) == len(
        variables
    )
    assert len(t1.functions) == len(t2.functions) and len(t1.functions) == len(
        functions
    )
    for v, f in zip(t1.variables, t1.functions):
        assert (
            (v in t2.variables)
            and (v in variables)
            and (f in t2.functions)
            and (f in functions)
        )

# Tests fot Tensor Class to support scalars and sympy expression type scalars

@pytest.mark.parametrize("scalar", [11.89, y * z + 5])
def test_tensor_scalar(scalar):
    scalar_tensor = Tensor(scalar)
    assert scalar_tensor.tensor().rank() == 0

# Tests for lambdify

@xfail
def test_lambdify_on_schwarzschild_metric_without_args():
    sch = schwarzschild_metric()
    
    # values of t, r, theta, phi, a, c
    
    vals = (0.0, 3.0, np.pi / 2, np.pi / 3, 2, 1)
    f = sch.tensor_lambdify()[1]
    
    # print(sch.variables)
    
    result_arr = np.array(f(*vals))
    cmp_arr = np.zeros((4, 4), dtype=float)
    cmp_arr[0, 0], cmp_arr[1, 1], cmp_arr[2, 2], cmp_arr[3, 3] = (
        1 - (vals[4] / vals[1]),
        -1 / ((1 - (vals[4] / vals[1])) * (vals[5] ** 2)),
        -1 * (vals[1] ** 2) / (vals[5] ** 2),
        -1 * (vals[1] ** 2) * (np.sin(vals[2]) ** 2) / (vals[5] ** 2),
    )
    
    assert_allclose(cmp_arr, result_arr, atol=1e-7, rtol=0.0)

@xfail
def test_lambdify_with_args():
    x, y = symbols("x y")
    T = BaseRelativityTensor([x + y, x], (x, y), config="l")
    args, f = T.tensor_lambdify(y, x)
    arr = np.array(f(2, 1))
    cmp_arr = np.array([3, 1])
    assert_allclose(arr, cmp_arr, rtol=0.0, atol=1e-7)
    for e1, e2 in zip(args, (y, x)):
        assert simplify(e1 - e2) == 0
