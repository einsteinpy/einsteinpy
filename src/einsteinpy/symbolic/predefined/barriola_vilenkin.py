from sympy import diag, sin, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def BarriolaVilekin(c=constants.c, k=symbols("k")):
    """
    Barriola-Vilekin monopol metric
    Phys. Rev. Lett. 63, 341
    Manuel Barriola and Alexander Vilenkin
    Published 24 July 1989

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.
    k : ~sympy.core.basic.Basic or int or float
        The scaling factor responsible for the deficit/surplus angle
        Defaults to ``k``.
    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    # define the metric
    metric = diag(
        -1, 1 / (c**2), ((k * r) ** 2) / (c**2), ((k * r * sin(th)) ** 2) / (c**2)
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="BarriolaVilekinMetric")
