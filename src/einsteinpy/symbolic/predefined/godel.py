from sympy import Rational, diag, exp, sqrt, symbols

from einsteinpy.symbolic.metric import MetricTensor


def Godel():
    """
    Godel metric
    Rev. Mod. Phys., v21, p447, (1949)
    Stephani (10.25) 122
    """
    coords = symbols("t x y z")
    om = symbols("omega")
    t, x, y, z = coords
    # define the metric
    metric = diag(-1, 1, -Rational(1, 2) * exp(2 * sqrt(2) * om * x), 1)
    metric[0, 2] = metric[2, 0] = -exp(sqrt(2) * om * x)
    metric = metric.tolist()
    return MetricTensor(metric, coords, "ll", name="GodelMetric")
