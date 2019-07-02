from sympy import Matrix, diag, symbols, sin, simplify
from einsteinpy.symbolic.tensor import *
from einsteinpy.symbolic.metric import *

_t, _r, _th, _ph = symbols("t r theta phi", real=True)
_coords = [_t, _r, _th, _ph]
_schw = diag(1 - 1 / _r, -1 / (1 - 1 / _r), -_r ** 2, -_r ** 2 * sin(_th) ** 2)
_E, _p1, _p2, _p3 = symbols("E p_1:4", positive=True)
_momentum = [_E, _p1, _p2, _p3]


def test_DiffOperator():
    pass


def test_PartialDerivative():
    pass
