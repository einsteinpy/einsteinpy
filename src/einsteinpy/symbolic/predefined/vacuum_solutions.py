from sympy import diag, sin, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def Schwarzschild(c=constants.c, sch=symbols("r_s")):
    """
    Schwarzschild exterior metric in curvature coordinates
    Schwarzschild, Sitz. Preuss. Akad. Wiss., p189, (1916)
    Stephani (13.19) p157

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.
    sch : ~sympy.core.basic.Basic or int or float
        Any value to assign to Schwarzschild Radius of the central object. 
        Defaults to 'r_s'.
    
    """
    coords = symbols("t r theta phi")
    t, r, theta, phi = coords
    val1, c2 = 1 - sch / r, c ** 2
    metric = diag(
        val1, -1 / (val1 * c2), -1 * (r ** 2) / c2, -1 * ((r * sin(theta)) ** 2) / c2
    ).tolist()
    return MetricTensor(metric, coords, "ll")
