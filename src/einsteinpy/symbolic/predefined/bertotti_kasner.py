from sympy import diag, exp, sin, sqrt, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def BertottiKasner(c=constants.c, k=symbols("k"), lambd=symbols("l")):
    """
    Birkhoff's theorem with Λ-term and Bertotti-Kasner space
    Phys. Lett. A, 245:363-365, 1998
    W. Rindler

    Parameters
    ----------
    c : ~sympy.core.basic.Basic or int or float
        Any value to assign to speed of light. Defaults to 'c'.
    lambd : ~sympy.core.basic.Basic or int or float
        The cosmological constant, note it must be postive.
        Defaults to ``l``.
    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    # define the metric
    metric = diag(
        -1,
        exp(2 * sqrt(lambd) * c * t) / (c**2),
        1 / (lambd * (c**2)),
        (sin(th) ** 2) / (lambd * (c**2)),
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="BertottiKasnerMetric")
