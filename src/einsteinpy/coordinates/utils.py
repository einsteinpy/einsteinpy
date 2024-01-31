import numpy as np

from einsteinpy import constant
from einsteinpy.ijit import jit

_c = constant.c.value


def calculate_coordinates(r, th, p, alpha):
    xa = np.sqrt(r**2 + alpha**2)
    sin_norm = xa * np.sin(th)
    e0 = sin_norm * np.cos(p)
    e1 = sin_norm * np.sin(p)
    e2 = r * np.cos(th)
    return xa, e0, e1, e2


def cartesian_to_spherical_fast(e0, e1, e2, e3, u0=None, u1=None, u2=None):
    vel = [u0, u1, u2]
    if all(u is not None for u in vel) and any(u != 0 for u in vel):
        return cartesian_to_spherical(e0, e1, e2, e3, u0, u1, u2)
    return cartesian_to_spherical_novel(e0, e1, e2, e3)


@jit
def cartesian_to_spherical(e0, e1, e2, e3, u0, u1, u2):
    """
    Utility function (jitted) to convert cartesian to spherical.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e4 = np.hypot(e1, e2)
    e5 = np.hypot(e4, e3)
    e6 = np.arctan2(e4, e3)
    e7 = np.arctan2(e2, e1)
    n1 = e1**2 + e2**2
    n2 = n1 + e3**2
    u3 = (e1 * u0 + e2 * u1 + e3 * u2) / np.sqrt(n2)
    u4 = (e3 * (e1 * u0 + e2 * u1) - n1 * u2) / (n2 * np.sqrt(n1))
    u5 = -1 * (u0 * e2 - e1 * u1) / n1

    return e0, e5, e6, e7, u3, u4, u5


@jit
def cartesian_to_spherical_novel(e0, e1, e2, e3):
    """
    Utility function (jitted) to convert cartesian to spherical.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e4 = np.hypot(e1, e2)
    e5 = np.hypot(e4, e3)
    e6 = np.arctan2(e4, e3)
    e7 = np.arctan2(e2, e1)

    return e0, e5, e6, e7


def cartesian_to_bl_fast(e0, e1, e2, e3, e4, u0=None, u1=None, u2=None):
    vel = [u0, u1, u2]
    if all(u is not None for u in vel) and any(u != 0 for u in vel):
        return cartesian_to_bl(e0, e1, e2, e3, e4, u0, u1, u2)
    return cartesian_to_bl_novel(e0, e1, e2, e3, e4)


@jit
def cartesian_to_bl(e0, e1, e2, e3, e4, u0, u1, u2):
    """
    Utility function (jitted) to convert cartesian to boyer lindquist.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e5 = (e1**2 + e2**2 + e3**2) - (e4**2)
    e6 = np.sqrt(0.5 * (e5 + np.sqrt((e5**2) + (4 * (e4**2) * (e3**2)))))
    e7 = np.arccos(e3 / e6)
    e8 = np.arctan2(e2, e1)
    e9 = 2 * (e1 * u0 + e2 * u1 + e3 * u2)
    u3 = (1 / (2 * e6)) * (
        (e9 / 2)
        + (
            (e5 * e9 + 4 * (e4**2) * e3 * u2)
            / (2 * np.sqrt((e5**2) + (4 * (e4**2) * (e3**2))))
        )
    )
    u4 = (-1 / np.sqrt(1 - np.square(e3 / e6))) * ((u2 * e6 - u3 * e3) / (e6**2))
    u5 = (1 / (1 + np.square(e2 / e1))) * ((u1 * e1 - u0 * e2) / (e1**2))

    return e0, e6, e7, e8, u3, u4, u5


@jit
def cartesian_to_bl_novel(e0, e1, e2, e3, e4):
    """
    Utility function (jitted) to convert cartesian to boyer lindquist.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e5 = (e1**2 + e2**2 + e3**2) - (e4**2)
    e6 = np.sqrt(0.5 * (e5 + np.sqrt((e5**2) + (4 * (e4**2) * (e3**2)))))
    e7 = np.arccos(e3 / e6)
    e8 = np.arctan2(e2, e1)

    return e0, e6, e7, e8


def spherical_to_cartesian_fast(e0, e1, e2, e3, u0=None, u1=None, u2=None):
    vel = [u0, u1, u2]
    if all(u is not None for u in vel) and any(u != 0 for u in vel):
        return spherical_to_cartesian(e0, e1, e2, e3, u0, u1, u2)
    return spherical_to_cartesian_novel(e0, e1, e2, e3)


@jit
def spherical_to_cartesian(e0, e1, e2, e3, u0, u1, u2):
    """
    Utility function (jitted) to convert spherical to cartesian.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e4 = e1 * np.cos(e3) * np.sin(e2)
    e5 = e1 * np.sin(e3) * np.sin(e2)
    e6 = e1 * np.cos(e2)
    u3 = (
        np.sin(e2) * np.cos(e3) * u0
        - e1 * np.sin(e2) * np.sin(e3) * u2
        + e1 * np.cos(e2) * np.cos(e3) * u1
    )
    u4 = (
        np.sin(e2) * np.sin(e3) * u0
        + e1 * np.cos(e2) * np.sin(e3) * u1
        + e1 * np.sin(e2) * np.cos(e3) * u2
    )
    u5 = np.cos(e2) * u0 - e1 * np.sin(e2) * u1

    return e0, e4, e5, e6, u3, u4, u5


@jit
def spherical_to_cartesian_novel(e0, e1, e2, e3):
    """x = r * np.cos(p) * np.sin(th)
    y = r * np.sin(p) * np.sin(th)
    z = r * np.cos(
    Utility function (jitted) to convert spherical to cartesian.
    This function should eventually result in Coordinate Transformation Graph!

    """
    e4 = e1 * np.cos(e3) * np.sin(e2)
    e5 = e1 * np.sin(e3) * np.sin(e2)
    e6 = e1 * np.cos(e2)

    return e0, e4, e5, e6


def bl_to_cartesian_fast(e0, e1, e2, e3, e4, u0=None, u1=None, u2=None):
    vel = [u0, u1, u2]
    if all(u is not None for u in vel) and any(u != 0 for u in vel):
        return bl_to_cartesian(e0, e1, e2, e3, e4, u0, u1, u2)
    return bl_to_cartesian_novel(e0, e1, e2, e3, e4)


@jit
def bl_to_cartesian(e0, e1, e2, e3, e4, u0, u1, u2):
    """
    Utility function (jitted) to convert bl to cartesian.
    This function should eventually result in Coordinate Transformation Graph!

    """
    xa, e5, e6, e7 = calculate_coordinates(e1, e2, e3, e4)

    u3 = (
        (e1 * u0 * np.sin(e2) * np.cos(e3) / xa)
        + (xa * np.cos(e2) * np.cos(e3) * u1)
        - (xa * np.sin(e2) * np.sin(e3) * u2)
    )
    u4 = (
        (e1 * u0 * np.sin(e2) * np.sin(e3) / xa)
        + (xa * np.cos(e2) * np.sin(e3) * u1)
        + (xa * np.sin(e2) * np.cos(e3) * u2)
    )
    u5 = (u0 * np.cos(e2)) - (e1 * np.sin(e2) * u1)

    return e0, e5, e6, e7, u3, u4, u5


@jit
def bl_to_cartesian_novel(e0, e1, e2, e3, e4):
    """
    Utility function (jitted) to convert bl to cartesian.
    This function should eventually result in Coordinate Transformation Graph!

    """
    _, e5, e6, e7 = calculate_coordinates(e1, e2, e3, e4)
    return e0, e5, e6, e7


def lorentz_factor(v1, v2, v3):
    """
    Returns the Lorentz Factor, ``gamma``

    Parameters
    ----------
    v1 : float
        First component of 3-Velocity
    v2 : float
        Second component of 3-Velocity
    v3 : float
        Third component of 3-Velocity

    Returns
    -------
    gamma : float
        Lorentz Factor

    """
    v_vec = np.array([v1, v2, v3])
    v_norm2 = v_vec.dot(v_vec)
    gamma = 1 / np.sqrt(1 - v_norm2 / _c**2)

    return gamma


@jit
def v0(g_cov_mat, v1, v2, v3):
    """
    Utility function to return Timelike component (v0) of 4-Velocity
    Assumes a (+, -, -, -) Metric Signature

    Parameters
    ----------
    g_cov_mat : ~numpy.ndarray
        Matrix, containing Covariant Metric \
        Tensor values, in same coordinates as ``v_vec``
        Numpy array of shape (4,4)
    v1 : float
        First component of 3-Velocity
    v2 : float
        Second component of 3-Velocity
    v3 : float
        Third component of 3-Velocity
    Returns
    -------
    float
        Timelike component of 4-Velocity

    """
    g = g_cov_mat
    # Factor to add to coefficient, C
    fac = -1 * _c**2
    # Defining coefficients for quadratic equation
    A = g[0, 0]
    B = 2 * (g[0, 1] * v1 + g[0, 2] * v2 + g[0, 3] * v3)
    C = (
        (g[1, 1] * v1**2 + g[2, 2] * v2**2 + g[3, 3] * v3**2)
        + 2 * v1 * (g[1, 2] * v2 + g[1, 3] * v3)
        + 2 * v2 * g[2, 3] * v3
        + fac
    )
    D = (B**2) - (4 * A * C)

    v_t = (-B + np.sqrt(D)) / (2 * A)

    return v_t
