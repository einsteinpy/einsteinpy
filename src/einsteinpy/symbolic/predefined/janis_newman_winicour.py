from sympy import diag, sin, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def JanisNewmanWinicour(
    c=constants.c, G=constants.G, gam=symbols("gam"), M=symbols("M")
):
    """
    Reality of the Schwarzschild singularity.
    Phys. Rev. Lett., 20:878-880, 1968.
    A. I. Janis, E. T. Newman, and J. Winicour.

    Parameters
    ----------
    M : ~sympy.core.basic.Basic or int or float
        Mass parameter, this is used for defining the schwarzschild metric. Defaults to ``M``.
    gam : ~sympy.core.basic.Basic or int or float
        Parameter for scaling  Schwarzschild radius, for gamma=1 this will return the  Schwarzschild metric
        Defaults to ``gam``.
    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    # Helper functions
    r_s = (2 * G * M) / (c**2)
    alpha = 1 - (r_s / (gam * r))

    # define the metric
    metric = diag(
        -1 * (alpha**gam),
        (alpha**-gam) / (c**2),
        (r**2) * (alpha ** (-gam + 1)),
        (r**2) * (alpha ** (-gam + 1)) * (sin(th) ** 2),
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="JanisNewmanWinicourMetric")
