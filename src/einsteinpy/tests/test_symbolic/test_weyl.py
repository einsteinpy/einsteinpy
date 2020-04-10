import numpy as np
import sympy
from sympy import cos, sin, sinh

from einsteinpy.symbolic import (
    MetricTensor,
    RicciScalar,
    RicciTensor,
    RiemannCurvatureTensor,
    WeylTensor,
    simplify_sympy_array,
)


def spherical_metric():
    symbolstr = "r theta phi xi"
    syms = sympy.symbols(symbolstr)
    list2d = (np.zeros(shape=(4, 4), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sympy.sin(syms[1]) ** 2)
    list2d[3][3] = (syms[0] ** 2) * (sympy.cos(syms[1]) ** 2)
    metric = MetricTensor(list2d, syms)
    return metric


def anti_de_sitter_metric():
    coords = sympy.symbols("t chi theta phi")
    t, ch, th, ph = coords
    m = sympy.diag(
        -1,
        cos(t) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2 * sin(th) ** 2,
    ).tolist()
    metric = MetricTensor(m, coords)
    return metric


def test_weyl_TypeError():
    testarr = np.ones((4, 4, 4, 4), dtype=int).tolist()
    syms = 0
    try:
        obj = WeylTensor(testarr, syms)
        assert False
    except TypeError:
        assert True


def test_weyl_ValueError_wrong_config_length():
    testarr = np.ones((4, 4, 4, 4), dtype=int).tolist()
    syms = sympy.symbols("x y z w")
    try:
        obj = WeylTensor(testarr, syms, config="uuu", parent_metric=None)
        assert False
    except ValueError:
        assert True


def test_weyl_dim2():
    list_metric = np.zeros((2, 2), dtype=int).tolist()
    symbolstr = "t r"
    syms = sympy.symbols(symbolstr)
    metric = MetricTensor(list_metric, syms)
    try:
        obj = WeylTensor.from_metric(metric)
        assert False
    except ValueError:
        assert True


def test_weyl_dim3():
    list_metric = np.ones((3, 3), dtype=int).tolist()
    symbolstr = "t r theta"
    syms = sympy.symbols(symbolstr)
    metric = MetricTensor(list_metric, syms)
    obj = WeylTensor.from_metric(metric).arr
    assert obj == sympy.Array(np.zeros((3, 3, 3, 3), dtype=int))


def test_weyl_conformal_rescaling():
    # https://en.wikipedia.org/wiki/Weyl_tensor#Conformal_rescaling
    a = sympy.symbols("a")
    mw1 = anti_de_sitter_metric()
    mw2 = MetricTensor(a * mw1.tensor(), mw1.symbols(), mw1.config)
    w1 = WeylTensor.from_metric(mw1).change_config("ulll")
    w2 = WeylTensor.from_metric(mw2).change_config("ulll")
    cmp_arr = sympy.Array(np.zeros(shape=w1.tensor().shape, dtype=int))
    assert simplify_sympy_array(w1.tensor() - w2.tensor()) == cmp_arr
    assert w1.syms == w1.symbols()
    assert w1.parent_metric == w1._parent_metric


def test_weyl_contraction_1st_3rd_indices_zero():
    mw1 = anti_de_sitter_metric()
    w1 = WeylTensor.from_metric(mw1)
    t1 = simplify_sympy_array(
        sympy.tensorcontraction(w1.change_config("ulll").tensor(), (0, 2))
    )
    t0 = sympy.Array(np.zeros(shape=t1.shape, dtype=int))
    assert t1 == t0
