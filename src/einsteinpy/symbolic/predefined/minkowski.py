from sympy import diag, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def Minkowski(c=constants.c):
    """
    Minkowski(flat) space-time. Space-time without any curvature or matter.

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.

    
    """
    coords = symbols("t x y z")
    metric = diag(-1, 1 / (c ** 2), 1 / (c ** 2), 1 / (c ** 2)).tolist()
    return MetricTensor(metric, coords, "ll", name="MinkowskiMetric")
