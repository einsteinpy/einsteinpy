import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.ijit import jit

_c = constant.c.value


def cartesian_to_spherical_fast(
    x, y, z, v_x=None, v_y=None, v_z=None, velocities_provided=False
):
    if velocities_provided:
        return cartesian_to_spherical(x, y, z, v_x, v_y, v_z)
    return cartesian_to_spherical_novel(x, y, z)


@jit
def cartesian_to_spherical(x, y, z, v_x, v_y, v_z):
    """
    Utility function (jitted) to convert cartesian to spherical.
    This function should eventually result in Coordinate Transformation Graph!
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arctan2(hxy, z)
    phi = np.arctan2(y, x)
    n1 = x ** 2 + y ** 2
    n2 = n1 + z ** 2
    v_r = (x * v_x + y * v_y + z * v_z) / np.sqrt(n2)
    v_t = (z * (x * v_x + y * v_y) - n1 * v_z) / (n2 * np.sqrt(n1))
    v_p = -1 * (v_x * y - x * v_y) / n1

    return r, theta, phi, v_r, v_t, v_p


@jit
def cartesian_to_spherical_novel(x, y, z):
    """
    Utility function (jitted) to convert cartesian to spherical.
    This function should eventually result in Coordinate Transformation Graph!
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arctan2(hxy, z)
    phi = np.arctan2(y, x)

    return r, theta, phi


def cartesian_to_bl_fast(
    x, y, z, a, v_x=None, v_y=None, v_z=None, velocities_provided=False
):
    if velocities_provided:
        return cartesian_to_bl(x, y, z, a, v_x, v_y, v_z)
    return cartesian_to_bl_novel(x, y, z, a)


@jit
def cartesian_to_bl(x, y, z, a, v_x, v_y, v_z):
    """
    Utility function (jitted) to convert cartesian to boyer lindquist.
    This function should eventually result in Coordinate Transformation Graph!
    """
    w = (x ** 2 + y ** 2 + z ** 2) - (a ** 2)
    r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (z ** 2)))))
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    dw_dt = 2 * (x * v_x + y * v_y + z * v_z)
    v_r = (1 / (2 * r)) * (
        (dw_dt / 2)
        + (
            (w * dw_dt + 4 * (a ** 2) * z * v_z)
            / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (z ** 2))))
        )
    )
    v_t = (-1 / np.sqrt(1 - np.square(z / r))) * ((v_z * r - v_r * z) / (r ** 2))
    v_p = (1 / (1 + np.square(y / x))) * ((v_y * x - v_x * y) / (x ** 2))

    return r, theta, phi, v_r, v_t, v_p, a


@jit
def cartesian_to_bl_novel(x, y, z, a):
    """
    Utility function (jitted) to convert cartesian to boyer lindquist.
    This function should eventually result in Coordinate Transformation Graph!
    """
    w = (x ** 2 + y ** 2 + z ** 2) - (a ** 2)
    r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (z ** 2)))))
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return r, theta, phi, a


def spherical_to_cartesian_fast(
    r, t, p, v_r=None, v_t=None, v_p=None, velocities_provided=False
):
    if velocities_provided:
        return spherical_to_cartesian(r, t, p, v_r, v_t, v_p)
    return spherical_to_cartesian_novel(r, t, p)


@jit
def spherical_to_cartesian(r, t, p, v_r, v_t, v_p):
    """
    Utility function (jitted) to convert spherical to cartesian.
    This function should eventually result in Coordinate Transformation Graph!
    """
    x = r * np.cos(p) * np.sin(t)
    y = r * np.sin(p) * np.sin(t)
    z = r * np.cos(t)
    v_x = (
        np.sin(t) * np.cos(p) * v_r
        - r * np.sin(t) * np.sin(p) * v_p
        + r * np.cos(t) * np.cos(p) * v_t
    )
    v_y = (
        np.sin(t) * np.sin(p) * v_r
        + r * np.cos(t) * np.sin(p) * v_t
        + r * np.sin(t) * np.cos(p) * v_p
    )
    v_z = np.cos(t) * v_r - r * np.sin(t) * v_t
    return x, y, z, v_x, v_y, v_z


@jit
def spherical_to_cartesian_novel(r, t, p):
    """
    Utility function (jitted) to convert spherical to cartesian.
    This function should eventually result in Coordinate Transformation Graph!
    """
    x = r * np.cos(p) * np.sin(t)
    y = r * np.sin(p) * np.sin(t)
    z = r * np.cos(t)
    return x, y, z


def bl_to_cartesian_fast(
    r, t, p, a, v_r=None, v_t=None, v_p=None, velocities_provided=False
):
    if velocities_provided:
        return bl_to_cartesian(r, t, p, a, v_r, v_t, v_p)
    return bl_to_cartesian_novel(r, t, p, a)


@jit
def bl_to_cartesian(r, t, p, a, v_r, v_t, v_p):
    """
    Utility function (jitted) to convert bl to cartesian.
    This function should eventually result in Coordinate Transformation Graph!
    """
    xa = np.sqrt(r ** 2 + a ** 2)
    sin_norm = xa * np.sin(t)
    x = sin_norm * np.cos(p)
    y = sin_norm * np.sin(p)
    z = r * np.cos(t)
    v_x = (
        (r * v_r * np.sin(t) * np.cos(p) / xa)
        + (xa * np.cos(t) * np.cos(p) * v_t)
        - (xa * np.sin(t) * np.sin(p) * v_p)
    )
    v_y = (
        (r * v_r * np.sin(t) * np.sin(p) / xa)
        + (xa * np.cos(t) * np.sin(p) * v_t)
        + (xa * np.sin(t) * np.cos(p) * v_p)
    )
    v_z = (v_r * np.cos(t)) - (r * np.sin(t) * v_t)
    return x, y, z, v_x, v_y, v_z


@jit
def bl_to_cartesian_novel(r, t, p, a):
    """
    Utility function (jitted) to convert bl to cartesian.
    This function should eventually result in Coordinate Transformation Graph!
    """
    xa = np.sqrt(r ** 2 + a ** 2)
    sin_norm = xa * np.sin(t)
    x = sin_norm * np.cos(p)
    y = sin_norm * np.sin(p)
    z = r * np.cos(t)
    return x, y, z


@jit
def lorentz_factor(v_vec):
    """
    Returns the Lorentz Factor, ``gamma``

    Parameters
    ----------
    v_vec : ~numpy.array
        Velocity 3-Vector

    Returns
    -------
    gamma : float
        Lorentz Factor

    """
    # Square of 3-Vector length
    v_norm2 = v_vec.dot(v_vec)
    gamma = 1 / np.sqrt(1 - v_norm2 / _c ** 2)

    return gamma


def v_t(g_cov_mat, v_vec, time_like=True):
    """
    Utility function to return Timelike component of 4-Velocity
    Assumes a (+, -, -, -) Metric Signature

    Parameters
    ----------
    g_cov_mat : ~numpy.array
        Matrix, containing Covariant Metric \
        Tensor values, in same coordinates as ``v_vec``
        Numpy array of shape (4,4)
    v_vec : ~numpy.array
        Velocity 3-Vector
    time_like : bool
        To determine, if the 4-Velocity is for a Time-like or \
        a Null-like Geodesic
        Defaults to ``True``

    Returns
    -------
    v_t : float
        Timelike component of 4-Velocity

    """
    u1, u2, u3 = v_vec
    g = g_cov_mat
    # Factor to add to ceofficient, C
    fac = -1 if time_like else 0
    # Defining coefficients for quadratic equation
    A = g[0, 0]
    B = 2 * (g[0, 1] * u1 + g[0, 2] * u2 + g[0, 3] * u3)
    C = (
        (g[1, 1] * u1 ** 2 + g[2, 2] * u2 ** 2 + g[3, 3] * u3 ** 2) + \
        2 * u1 * (g[1, 2] * u2 + g[1, 3] * u3) + \
        2 * u2 * g[2, 3] * u3 + \
        fac
    )
    D = (B ** 2) - (4 * A * C)

    v_t = (-B + np.sqrt(D)) / (2 * A)

    return v_t * u.one


def four_position(t, x_vec):
    """
    Utility function to return 4-Position

    Parameters
    ----------
    t : float
        Coordinate Time
    x_vec : ~numpy.array
        Position 3-Vector

    Returns
    -------
    x_4vec : ~numpy.array
        Position 4-Vector

    """
    x_4vec = np.append([_c * t], x_vec)

    return x_4vec


def four_velocity(g_cov_mat, v_vec, time_like):
    """
    Utility function to return 4-Velocity

    Parameters
    ----------
    g_cov_mat : ~numpy.array
        Matrix, containing Covariant Metric \
        Tensor values, in same coordinates as ``v_vec``
        Numpy array of shape (4,4)
    v_vec : ~numpy.array
        Velocity 3-Vector
    time_like : bool
        To determine, if the 4-Velocity is for a Time-like or \
        a Null-like Geodesic
        Defaults to ``True``

    Returns
    -------
    v_4vec : ~numpy.array
        Velocity 4-Vector

    """
    v_4vec = np.append(v_t(g_cov_mat, v_vec, time_like), v_vec)

    return v_4vec


def stacked_vec(g_cov_mat, t, x_vec, v_vec, time_like):
    """
    Packages 4-Position and 4-Velocity into a Length-8 vector

    Parameters
    ----------
    g_cov_mat : ~numpy.array
        Matrix, containing Covariant Metric \
        Tensor values, in same coordinates as \
        ``x_vec`` or ``v_vec``
        Numpy array of shape (4,4)
    t : float
        Coordinate Time
    x_vec : ~numpy.array
        Position 3-Vector
    v_vec : ~numpy.array
        Velocity 3-Vector
    time_like : bool
        To determine, if the 4-Velocity is for a Time-like or \
        a Null-like Geodesic
        Defaults to ``True``

    Returns
    -------
    stacked_vec : ~numpy.array
        Length-8 Vector of form [x0, x1, x2, x3, v0, v1, v2, v3]

    """
    x_4vec = four_position(t, x_vec)
    v_4vec = four_velocity(g_cov_mat, v_vec, time_like)

    stacked_vec = np.hstack((x_4vec, v_4vec))

    return stacked_vec
