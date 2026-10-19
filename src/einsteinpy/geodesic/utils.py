"""
Utilities for Geodesic Module

Unit System: M-Units => :math:`c = G = M = k_e = 1`
Metric Signature => :math:`(-, +, +, +)`

"""
import numpy as np

from einsteinpy.utils.dual import DualNumber


def _P(g, g_prms, q, p, time_like=True):
    """
    Utility function to compute 4-Momentum of the test particle

    Parameters
    ----------
    g : callable
        Metric Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    q : array_like
        Initial 4-Position
    p : array_like
        Initial 3-Momentum
    time_like: bool, optional
        Determines type of Geodesic
        ``True`` for Time-like geodesics
        ``False`` for Null-like geodesics
        Defaults to ``True``

    Returns
    -------
    P: numpy.ndarray
        4-Momentum

    """
    guu = g(q, *g_prms)
    P = np.array([0.0, *p])

    A = guu[0, 0]
    B = 2 * guu[0, 3] * P[3]
    C = (
        guu[1, 1] * P[1] * P[1]
        + guu[2, 2] * P[2] * P[2]
        + guu[3, 3] * P[3] * P[3]
        + int(time_like)
    )

    P[0] = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)

    return P


def sigma(r, theta, a):
    """
    Returns the value of :math:`r^2 + a^2 * \\cos^2(\\theta)`
    Uses natural units, with :math:`c = G = M = k_e = 1`

    Parameters
    ----------
    r : float
        r-component of 4-Position
    theta : float
        theta-component of 4-Position
    a : float
        Spin Parameter
        :math:`0 \\le a \\le 1`

    Returns
    -------
    float
        The value of :math:`r^2 + a^2 * \\cos^2(\\theta)`

    """
    sigma = (r**2) + ((a * np.cos(theta)) ** 2)

    return sigma


def delta(r, a, Q=0):
    """
    Returns the value of :math:`r^2 - r_s r + a^2 + r_Q^2`
    Uses natural units, with :math:`c = G = M = k_e = 1`

    Parameters
    ----------
    r : float
        r-component of 4-Position
    a : float
        Spin Parameter
        :math:`0 \\le a \\le 1`
    Q : float
        Charge on gravitating body
        Defaults to ``0``

    Returns
    -------
    float
        The value of :math:`r^2 - r_s r + a^2 + r_Q^2`

    """
    delta = (r**2) - (2 * r) + (a**2) + Q**2

    return delta


def _sch(x_vec, *params):
    """
    Contravariant Schwarzschild Metric in Spherical Polar coordinates
    Uses natural units, with :math:`c = G = M = k_e = 1`

    Parameters
    ----------
    x_vec : array_like
        4-Position

    Other Parameters
    ----------------
    params : array_like
        Tuple of parameters to pass to the metric

    Returns
    -------
    numpy.ndarray
        Contravariant Schwarzschild Metric Tensor

    """
    r, th = x_vec[1], x_vec[2]

    g = np.zeros(shape=(4, 4), dtype=DualNumber)

    tmp = 1.0 - (2 / r)
    g[0, 0] = -1 / tmp
    g[1, 1] = tmp
    g[2, 2] = 1 / (r**2)
    g[3, 3] = 1 / ((r * np.sin(th)) ** 2)

    return g


def _kerr(x_vec, *params):
    """
    Contravariant Kerr Metric in Boyer-Lindquist coordinates
    Uses natural units, with :math:`c = G = M = k_e = 1`

    Parameters
    ----------
    x_vec : array_like
        4-Position

    Other Parameters
    ----------------
    params : array_like
        Tuple of parameters to pass to the metric
        Should contain Spin Parameter, ``a``

    Returns
    -------
    numpy.ndarray
        Contravariant Kerr Metric Tensor

    """
    a = params[0]

    r, th = x_vec[1], x_vec[2]
    sg, dl = sigma(r, th, a), delta(r, a)

    g = np.zeros(shape=(4, 4), dtype=DualNumber)

    g[0, 0] = -(r**2 + a**2 + (2 * r * (a * np.sin(th)) ** 2) / sg) / dl
    g[1, 1] = dl / sg
    g[2, 2] = 1 / sg
    g[3, 3] = (1 / (dl * np.sin(th) ** 2)) * (1 - 2 * r / sg)
    g[0, 3] = g[3, 0] = -(2 * r * a) / (sg * dl)

    return g


def _kerrnewman(x_vec, *params):
    """
    Contravariant Kerr-Newman Metric in Boyer-Lindquist coordinates
    Uses natural units, with :math:`c = G = M = k_e = 1`

    Parameters
    ----------
    x_vec : array_like
        4-Position

    Other Parameters
    ----------------
    params : array_like
        Tuple of parameters to pass to the metric
        Should contain Spin, ``a``, and Charge, ``Q``

    Returns
    -------
    numpy.ndarray
        Contravariant Kerr-Newman Metric Tensor

    """
    a, Q = params[0], params[1]

    r, th = x_vec[1], x_vec[2]
    sg, dl = sigma(r, th, a), delta(r, a, Q)
    a2 = a**2
    r2 = r**2
    sint2 = np.sin(th) ** 2
    csct2 = 1 / sint2
    csct4 = 1 / sint2**2

    g = np.zeros(shape=(4, 4), dtype=DualNumber)

    denom = dl * (a2 + 2 * r2 + a2 * np.cos(2 * th)) ** 2
    g[0, 0] = -(4 * sg * ((a2 + r2) ** 2 - a2 * dl * sint2) / denom)
    g[1, 1] = dl / sg
    g[2, 2] = 1 / sg
    g[3, 3] = (sg * csct4 * (-a2 + dl * csct2)) / (dl * (a2 - (a2 + r2) * csct2) ** 2)
    g[0, 3] = g[3, 0] = -(4 * a * (a2 - dl + r2) * sg / denom)

    return g
