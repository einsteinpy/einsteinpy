import numpy as np
from sympy import Matrix, diff, eye, sin, symbols, zeros


def christoffels(list2d, syms):
    """
    Function to calculate christoffel symbols of a given metric

    Parameters
    ----------
    list2d : list
        2d list (Matrix) representing metric, containing ~sympy expressions
    syms : list
        1d list containing representaion of [x0,x1,x2...] in ~sympy expressions

    Returns
    -------
    a : list
        3d list of ~sympy expressions containing christoffel symbols

    """
    dims = len(syms)
    christlist = (np.zeros(shape=(dims, dims, dims), dtype=int)).tolist()
    mat = Matrix(list2d)
    mat_inv = mat.inv()
    momocow = [i for i in range(dims * dims * dims)]
    for t in momocow:
        k = t % dims
        j = (int(t / dims)) % (dims)
        i = (int(t / (dims * dims))) % (dims)
        temp = 0
        for n in range(dims):
            temp += (mat_inv[i, n] / 2) * (
                diff(list2d[n][j], syms[k])
                + diff(list2d[n][k], syms[j])
                - diff(list2d[j][k], syms[n])
            )
        christlist[i][j][k] = temp
    return christlist


def schwarzschild_christoffels(symbolstr="t r theta phi"):
    """Returns the 3d list of christoffel symbols of Schwarzschild Metric."""
    list2d = [[0 for i in range(4)] for i in range(4)]
    syms = symbols(symbolstr)
    c, a = symbols("c a")
    list2d[0][0] = 1 - (a / syms[1])
    list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
    list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
    list2d[3][3] = -1 * (syms[1] ** 2) * (sin(syms[2]) ** 2) / (c ** 2)
    return christoffels(list2d, syms)
