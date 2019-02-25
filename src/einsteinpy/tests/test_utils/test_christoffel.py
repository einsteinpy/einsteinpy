import numpy as np
import pytest
import sympy

from einsteinpy.utils import christoffel


def test_schwarzschild_christoffels():
    symbolstr = "t r theta phi"
    syms = sympy.symbols(symbolstr)
    l = christoffel.schwarzschild_christoffels(symbolstr)
    G, M, c, a = sympy.symbols("G M c a")
    assert (
        sympy.simplify(l[2][3][3] - (-1 * sympy.cos(syms[2]) * sympy.sin(syms[2]))) == 0
    )
    assert sympy.simplify(l[3][3][1] - syms[1] / (syms[1] ** 2)) == 0
    assert (
        sympy.simplify(
            (l[1][1][1].subs({a: (2 * G * M / (c ** 2))}))
            - (G * M / (2 * G * M * syms[1] - c ** 2 * syms[1] ** 2))
        )
        == 0
    )


def test_riemann_curvature_tensor():
    symbolstr = "r theta phi"
    syms = sympy.symbols(symbolstr)
    list2d = (np.zeros(shape=(3, 3), dtype=int)).tolist()
    list2d[0][0] = 1
    list2d[1][1] = syms[0] ** 2
    list2d[2][2] = (syms[0] ** 2) * (sympy.sin(syms[1]) ** 2)
    T = christoffel.riemann_curvature_tensor(list2d, syms)
    Tnp = np.array(T)
    assert np.all(Tnp == 0)


def test_simplify_christoffels():
    a = sympy.symbols("a")
    b = a ** 2 - (1 / (1 / (a ** 3 / a))) - 1 + sympy.sin(a) ** 2 + sympy.cos(a) ** 2
    c = christoffel.simplify_christoffels([[[b]]], 1)
    assert c[0][0][0] == 0


def test_kerr_christoffels():
    a = christoffel.kerr_christoffels()
    assert a[0][0][0] == 0
