import numpy as np
from sympy import cos, simplify, symbols

from einsteinpy.symbolic import GenericVector, MetricTensor


def euclidean_space_metric():
    symbolstr = "e1 e2"  # let angle between e1 & e2 be theta
    syms = symbols(symbolstr)
    th = symbols("theta")
    list2d = np.zeros((2, 2), dtype=int).tolist()
    # defining the metric tensor when axis are not orthogonal
    list2d[0][0] = list2d[1][1] = 1
    list2d[1][0] = list2d[0][1] = cos(th)
    metric = MetricTensor(list2d, syms, config="ll")
    return metric


def test_GenericVector_change_config_theoretical_test():
    # https://en.wikipedia.org/wiki/Covariance_and_contravariance_of_vectors#Definition
    a, b, th = symbols("a b theta")
    metric = euclidean_space_metric()
    # defining a contravariant vector
    cnvec = GenericVector([a, b], metric.syms, config="u", parent_metric=metric)
    covec = cnvec.change_config("l")  # get contravariant vector
    assert simplify(covec.tensor()[0] - (a + b * cos(th))) == 0
    assert simplify(covec.tensor()[1] - (b + a * cos(th))) == 0


def test_GenericVector_check_ValueErrors():
    a, b = symbols("a b")
    syms = symbols("e1 e2")
    # input a tensor with wring rank
    try:
        arr = [[a, b], [b, 1]]
        v1 = GenericVector(arr, syms, "l")
        boolstore = False
    except ValueError:
        boolstore = True
    assert boolstore
    # input a wrong length config
    try:
        arr = [a, b]
        v2 = GenericVector(arr, syms, "uu")
        boolstore = False
    except ValueError:
        boolstore = True
    assert boolstore
