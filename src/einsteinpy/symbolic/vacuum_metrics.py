import sympy

from einsteinpy.symbolic.metric import MetricTensor


def SchwarzschildMetric(symbolstr="t r theta phi"):
    """
    Returns  Metric Tensor  of symbols of Schwarzschild Metric.

    Parameters
    ----------
    symbolstr : string
        symbols to be used to define schwarzschild space, defaults to 't r theta phi'

    Returns
    -------
    ~einsteinpy.symbolic.metric.MetricTensor
        Metric Tensor for Schwarzschild space-time
    """
    list2d = [[0 for i in range(4)] for i in range(4)]
    syms = sympy.symbols(symbolstr)
    c, a = sympy.symbols("c a")
    list2d[0][0] = 1 - (a / syms[1])
    list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
    list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
    list2d[3][3] = -1 * (syms[1] ** 2) * (sympy.sin(syms[2]) ** 2) / (c ** 2)

    return MetricTensor(list2d, syms)
