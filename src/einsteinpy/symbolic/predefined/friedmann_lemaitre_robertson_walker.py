from sympy import symbols, sin, diag, Function
from einsteinpy.symbolic.metric import MetricTensor

def FLRW():
    """
    Friedmann-Lemaître-Robertson-Walker (FLRW) metric
    in comoving spherical coordinates (t, r, theta, phi).

    Signature: (-, +, +, +).

    ds^2 = -dt^2 + a(t)^2 [ dr^2 / (1 - k r^2) + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2 ]

    Reference:
        Lachièze-Rey, M., & Luminet, J.-P. (1995).
        "Cosmic topology."
        Physics Reports, 254(3), 135–214.
        doi:10.1016/0370-1573(94)00085-H
    """
    # Coordinates and parameters
    t, r, th, ph = symbols("t r theta phi", real=True)
    k = symbols("k", real=True)              # curvature parameter (-1, 0, +1)
    a = Function("a")(t)                     # scale factor a(t)

    # Metric diagonal components
    g_tt   = -1
    g_rr   = a**2 / (1 - k*r**2)
    g_thth = a**2 * r**2
    g_phph = a**2 * r**2 * sin(th)**2

    metric = diag(g_tt, g_rr, g_thth, g_phph).tolist()
    coords = (t, r, th, ph)
    return MetricTensor(metric, coords, "ll", name="FLRWMetric")
