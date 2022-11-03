from sympy import diag, exp, sin, sqrt, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def Ernst(B=symbols("B"), M=symbols("M")):
    """
    Black holes in a magnetic universe.
    J. Math. Phys., 17:54-56, 1976.
    Frederick J. Ernst.

    Parameters
    ----------
    M : ~sympy.core.basic.Basic or int or float
        Mass of the black hole. Defaults to ``M``.
    B : ~sympy.core.basic.Basic or int or float
        The magnetic field strength
        Defaults to ``B``.
    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    # Helper functions
    lambd = 1 + ((B * r * sin(th)) ** 2)
    w = 1 - ((2 * M) / r)

    # define the metric
    metric = diag(
        -1 * (lambd**2) * w,
        (lambd**2) / w,
        ((r * lambd) ** 2),
        (((r * sin(th)) / lambd) ** 2),
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="ErnstMetric")
