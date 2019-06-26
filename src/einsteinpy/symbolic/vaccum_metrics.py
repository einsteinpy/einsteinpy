import sympy

from .metric import MetricTensor


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


def KerrMetric(symbolstr="t r theta phi"):
    """
    Returns Metric Tensor of symbols of Kerr metric(BL coordinates) in Plank units : G=1, c=1.
    Parameters
    ----------
    symbolstr : string
        symbols to be used to define kerr space in BL coordinates, defaults to 't r theta phi'
    Returns
    -------
    ~einsteinpy.symbolic.metric.MetricTensor
    Metric Tensor for Kerr space-time
    """
    list2d = [[0 for i in range(4)] for i in range(4)]
    syms = sympy.symbols(symbolstr)
    a, R = sympy.symbols("a R")
    A = syms[1] ** 2 - R * syms[1] + a ** 2
    sigma = syms[1] ** 2 + (a ** 2) * (sympy.cos(syms[2]) ** 2)
    list2d[0][0] = (R * syms[1] / sigma) - 1
    list2d[1][1] = sigma / A
    list2d[2][2] = sigma
    list2d[3][3] = (
        (sympy.sin(syms[2]) ** 2)
        * ((a ** 2 + syms[1] ** 2) ** 2 - (a ** 2) * (A * (sympy.sin(syms[2]) ** 2)))
    ) / sigma
    list2d[3][0] = -1 * (R * a * (syms[1])) * (sympy.sin(syms[2]) ** 2) / sigma
    list2d[0][3] = list2d[3][0]
    return MetricTensor(list2d, syms)
