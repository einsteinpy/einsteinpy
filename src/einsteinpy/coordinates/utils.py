import numpy as np

from einsteinpy.ijit import jit


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
