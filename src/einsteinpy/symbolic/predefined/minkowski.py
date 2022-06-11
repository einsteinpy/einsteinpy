from sympy import diag, sin, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def MinkowskiCartesian(c=constants.c):
    """
    Minkowski(flat) space-time in Cartesian coordinates. Space-time without any curvature or matter.

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.


    """
    coords = symbols("t x y z")
    metric = diag(-1, 1 / (c**2), 1 / (c**2), 1 / (c**2)).tolist()
    return MetricTensor(metric, coords, "ll", name="MinkowskiMetric")


Minkowski = MinkowskiCartesian


def MinkowskiPolar(c=constants.c):
    """
    Minkowski(flat) space-time in Polar coordinates. Space-time without any curvature or matter.

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.


    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    c2 = c**2
    metric = diag(-1, 1 / c2, (r**2) / c2, (r**2 * sin(th) ** 2) / c2).tolist()
    return MetricTensor(metric, coords, "ll", name="MinkowskiMetricPolar")
