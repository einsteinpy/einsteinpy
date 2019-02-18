import numpy as np
import pytest
from sympy import cos, simplify, sin, symbols

from einsteinpy.utils import christoffel


def test_schwarzschild_christoffels():
    symbolstr = "t r theta phi"
    syms = symbols(symbolstr)
    l = christoffel.schwarzschild_christoffels(symbolstr)
    G, M, c, a = symbols("G M c a")
    assert simplify(l[2][3][3] - (-1 * cos(syms[2]) * sin(syms[2]))) == 0
    assert simplify(l[3][3][1] - syms[1] / (syms[1] ** 2)) == 0
    assert (
        simplify(
            (l[1][1][1].subs({a: (2 * G * M / (c ** 2))}))
            - (G * M / (2 * G * M * syms[1] - c ** 2 * syms[1] ** 2))
        )
        == 0
    )


def test_riemann_curvature_tensor():
    symbolstr = "r theta phi"
    syms = symbols(symbolstr)
    list2d = (np.zeros(shape=(3, 3), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sin(syms[1]) ** 2)
    T = christoffel.riemann_curvature_tensor(list2d, syms)
    Tnp = np.array(T)
    assert np.all(Tnp == 0)
