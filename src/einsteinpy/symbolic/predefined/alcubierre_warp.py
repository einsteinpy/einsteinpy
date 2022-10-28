from sympy import diag, sqrt, symbols, tanh

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def AlcubierreWarp(
    x_s=symbols("x_s"), sigma=symbols("sigma"), R=symbols("R"), v=symbols("v")
):
    """
    Alcubierre Warp Drive Metric (G = c = 1) [1]_

    .. [1] Classical and Quantum Gravity,
    "The warp drive: hyper-fast travel within general relativity",
    Miguel Alcubierre, 1994, **11(5)**, pp. 73-77

    arXiv: `<https://arxiv.org/abs/gr-qc/0009013>`_)

    Journal: `<https://doi.org/10.1088%2F0264-9381%2F11%2F5%2F001>`_

    Parameters
    ----------
    x_s : ~sympy.core.basic.Basic or int or float
        Coordinate (a function of time, `t`)
    sigma : ~sympy.core.basic.Basic or int or float
        Arbitrary Parameter (See paper [1]_)
    R : ~sympy.core.basic.Basic or int or float
        Arbitrary Parameter (See paper [1]_)
    v : ~sympy.core.basic.Basic or int or float
        Coordinate Velocity

    Returns
    -------
    ~einsteinpy.symbolic.metric.MetricTensor
        Alcubierre Warp Drive Metric Tensor


    Alcubierre Warp Drive Metric (G = c = 1) [1]_

    """
    coords = symbols("t x y z")
    t, x, y, z = coords
    r_s = sqrt((x - x_s) ** 2 + y + z)
    f = (tanh(sigma * (r_s + R)) - tanh(sigma * (r_s - R))) / (2 * tanh(sigma * R))

    v2 = v**2
    f2 = f**2

    metric = diag(v2 * f2 - 1, 1, 1, 1).tolist()
    metric[0][1] = metric[1][0] = -v * f

    return MetricTensor(metric, coords, "ll", name="AlcubierreWarpMetric")
