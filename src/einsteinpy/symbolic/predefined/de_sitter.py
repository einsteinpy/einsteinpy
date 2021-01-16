from sympy import cos, cosh, diag, exp, sin, sinh, symbols

from einsteinpy.symbolic.metric import MetricTensor


def AntiDeSitter():
    """
    Anti-de Sitter space

    Hawking and Ellis (5.9) p131

    """
    coords = symbols("t chi theta phi")
    t, ch, th, ph = coords
    metric = diag(
        -1,
        cos(t) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2,
        cos(t) ** 2 * sinh(ch) ** 2 * sin(th) ** 2,
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="AntiDeSitterMetric")


def AntiDeSitterStatic():
    """
    Static form of Anti-de Sitter space

    Hawking and Ellis (5.9) p131

    """
    coords = symbols("t r theta phi")
    t, r, th, ph = coords
    metric = diag(-cosh(r) ** 2, 1, sinh(r) ** 2, sinh(r) ** 2 * sin(th) ** 2).tolist()
    return MetricTensor(metric, coords, "ll", name="AntiDeSitterStaticMetric")


def DeSitter():
    """
    de Sitter space

    Hawking and Ellis p125

    """
    coords = symbols("t x y z")
    t = coords[1]
    al = symbols("alpha")
    expr = exp(2 * t / al)
    metric = diag(-1, expr, expr, expr).tolist()
    return MetricTensor(metric, coords, "ll", name="DeSitterMetric")
