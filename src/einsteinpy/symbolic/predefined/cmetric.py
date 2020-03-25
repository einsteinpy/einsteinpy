from sympy import Function, diag, symbols

from einsteinpy.symbolic.metric import MetricTensor


def CMetric():
    """
    The C-metric
    Stephani (Table 16.2) p188
    """

    coords = symbols("t x y phi")
    x, y = coords[1], coords[2]
    f, h = Function("f")(x), Function("h")(y)
    metric = diag(
        -h / (x + y) ** 2,
        1 / ((x + y) ** 2 * f),
        1 / ((x + y) ** 2 * h),
        f / (x + y) ** 2,
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="CMetric")
