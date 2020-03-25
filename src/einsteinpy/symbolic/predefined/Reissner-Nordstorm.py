from sympy import cos, diag, pi, sin, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def Reissner_Nordstorm(
    c=constants.c,
    G=constants.G,
    eps_0=constants.eps_0,
    sch=symbols("r_s"),
    a=symbols("a"),
    Q=symbols("Q"),
):
    """
    The Reissner–Nordström metric in spherical coordinates

    A static solution to the Einstein–Maxwell field equations,
    which corresponds to the gravitational field of a charged,
    non-rotating, spherically symmetric body of mass M.

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to ``c``.
    G : ~sympy.core.basic.Basic or int or float
        Any value to assign to the Newton's (or gravitational) constant. Defaults to ``G``.
    eps_0 : ~sympy.core.basic.Basic or int or float
        Any value to assign to the electric constant or permittivity of free space. Defaults to ``eps_0``.
    sch : ~sympy.core.basic.Basic or int or float
        Any value to assign to Schwarzschild Radius of the central object.
        Defaults to ``r_s``.
    a : ~sympy.core.basic.Basic or int or float
        Spin factor of the heavy body. Usually, given by ``J/(Mc)``,
        where ``J`` is the angular momentum.
        Defaults to ``a``.
    Q:  ~sympy.core.basic.Basic or int or float
        Any value to assign to eletric charge of the central object.
        Defaults to ``Q``.

    """
    coords = symbols("t r theta phi")
    t, r, theta, phi = coords
    rQsq = ((Q ** 2) * G) / (4 * pi * eps_0 * (c ** 4))
    Arn = 1 - sch / r + rQsq / r ** 2

    metric = diag(
        (Arn * (c ** 2)), -(1 / Arn), -(r ** 2), -(r ** 2) * sin(theta) ** 2,
    ).tolist()
    return MetricTensor(metric, coords, "ll")
