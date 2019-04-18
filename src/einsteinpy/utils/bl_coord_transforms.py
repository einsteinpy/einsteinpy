import warnings

import astropy.units as u
import numpy as np

warnings.warn(
    "This module will be deprecated in EinsteinPy v0.3.0. Please switch to `einsteinpy.coordinates`.",
    PendingDeprecationWarning,
)


def CartesianToBL_pos(pos_vec, a):
    """
    Function to convert Cartesian to Boyer-Lindquits coordinates (position vector)

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x, y, z coordinates in SI units(m)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with r, theta, phi in (m, rad, rad)

    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    w = float(np.sum(np.square(pos_vec[:3])) - (a ** 2))
    r_vec[0] = np.sqrt(
        0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2))))
    )
    r_vec[1] = np.arccos(pos_vec[2] / r_vec[0])
    r_vec[2] = np.arctan2(pos_vec[1], pos_vec[0])
    return r_vec


def CartesianToBL_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocity in Cartesian coordinates to velocity in Boyer-Lindquist coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x, y, z coordinates in SI units(m)
    vel_vec : ~numpy.array
        3-length numpy array with vx, vy, vz in (m/s)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array having v_r, v_theta, v_phi in (m/s, rad/s, rad/s)

    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    w = float(np.sum(np.square(pos_vec[:3])) - (a ** 2))
    dw_dt = float(2 * np.sum(np.multiply(pos_vec, vel_vec)))
    r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2)))))
    v_vec[0] = (1 / (2 * r)) * (
        (dw_dt / 2)
        + (
            (w * dw_dt + 4 * (a ** 2) * pos_vec[2] * vel_vec[2])
            / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2))))
        )
    )
    v_vec[1] = (-1 / np.sqrt(1 - np.square(pos_vec[2] / r))) * (
        (vel_vec[2] * r - v_vec[0] * pos_vec[2]) / (r ** 2)
    )
    v_vec[2] = (1 / (1 + np.square(pos_vec[1] / pos_vec[0]))) * (
        (vel_vec[1] * pos_vec[0] - vel_vec[0] * pos_vec[1]) / (pos_vec[0] ** 2)
    )
    return v_vec


def BLToCartesian_pos(pos_vec, a):
    """
    Function to convert Boyer-Lindquist coordinates to Cartesian coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with x, y, z in m.

    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    tmp = np.sqrt(pos_vec[0] ** 2 + a ** 2) * np.sin(pos_vec[1])
    r_vec[0] = tmp * np.cos(pos_vec[2])
    r_vec[1] = tmp * np.sin(pos_vec[2])
    r_vec[2] = pos_vec[0] * np.cos(pos_vec[1])
    return r_vec


def BLToCartesian_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocities in Boyer-Lindquist coordinates to velocities in Cartesian coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    vel_vec : ~numpy.array
        3-length numpy array with V_r, V_theta, V_phi (m/s, rad/s, rad/s)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array having vx, vy, vz in SI units(m/s)

    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    tmp = np.sqrt(pos_vec[0] ** 2 + a ** 2)
    v_vec[0] = (
        (pos_vec[0] * vel_vec[0] * np.sin(pos_vec[1]) * np.cos(pos_vec[2]) / tmp)
        + (tmp * np.cos(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[1])
        - (tmp * np.sin(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[2])
    )
    v_vec[1] = (
        (pos_vec[0] * vel_vec[0] * np.sin(pos_vec[1]) * np.sin(pos_vec[2]) / tmp)
        + (tmp * np.cos(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[1])
        + (tmp * np.sin(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[2])
    )
    v_vec[2] = (vel_vec[0] * np.cos(pos_vec[1])) - (
        pos_vec[0] * np.sin(pos_vec[1]) * vel_vec[1]
    )
    return v_vec


def C2BL_units(pos_vec, vel_vec, a):
    """
    Function to convert Cartesian to Boyer-Lindquist Coordinates along with handling units

    Parameters
    ----------
    pos_vec : list
        list of 3 position coordinates along with appropriate units
        [x, y, z]
        (u.m, u.m, u.m)
    vel_vec : list
        list of 3 velocity coordinates along with appropriate units
        [vx, vy, vz]
        (u.m/u.s, u.m/u.s, u.m/u.s)
    a : float
        Any constant

    Returns
    -------
    tuple
        consisting of 2 lists
        ([r, theta, phi], [vr, vtheta, vphi]) in units
        ([u.m, u.rad, u.rad],[u.m/u.s, u.rad/u.s, u.rad/u.s])

    """
    units_list = [u.s, u.m, u.m, u.m, u.one, u.m / u.s, u.m / u.s, u.m / u.s]
    pos_vec_vals = [pos_vec[i].to(units_list[i + 1]).value for i in range(len(pos_vec))]
    vel_vec_vals = [vel_vec[i].to(units_list[i + 5]).value for i in range(len(vel_vec))]
    new_pos_units_list = [u.m, u.rad, u.rad]
    new_vel_units_list = [u.m / u.s, u.rad / u.s, u.rad / u.s]
    npos_vec = CartesianToBL_pos(np.array(pos_vec_vals), a).tolist()
    npos_vec = [e1 * e2 for e1, e2 in zip(npos_vec, new_pos_units_list)]
    nvel_vec = CartesianToBL_vel(np.array(pos_vec_vals), np.array(vel_vec_vals), a)
    nvel_vec = [e1 * e2 for e1, e2 in zip(nvel_vec, new_vel_units_list)]
    return (npos_vec, nvel_vec)


def BL2C_units(pos_vec, vel_vec, a):
    """
    Function to convert Boyer-Lindquist to Cartesian Coordinates along with handling units

    Parameters
    ----------
    pos_vec : list
        list of 3 position coordinates along with appropriate units
        [r, theta, phi]
        (u.m, u.rad, u.rad)
    vel_vec : list
        list of 3 velocity coordinates along with appropriate units
        [vr, vtheta, vphi]
        (u.m/u.s, u.rad/u.s, u.rad/u.s)
    a : float
        Any constant

    Returns
    -------
    tuple
        consisting of 2 lists
        ([x, y, z], [vx, vy, vz]) in units
        ([u.m, u.m, u.m],[u.m/u.s, u.m/u.s, u.m/u.s])

    """
    units_list = [u.s, u.m, u.rad, u.rad, u.one, u.m / u.s, u.rad / u.s, u.rad / u.s]
    pos_vec_vals = [pos_vec[i].to(units_list[i + 1]).value for i in range(len(pos_vec))]
    vel_vec_vals = [vel_vec[i].to(units_list[i + 5]).value for i in range(len(vel_vec))]
    new_pos_units_list = [u.m, u.m, u.m]
    new_vel_units_list = [u.m / u.s, u.m / u.s, u.m / u.s]
    npos_vec = BLToCartesian_pos(np.array(pos_vec_vals), a).tolist()
    npos_vec = [e1 * e2 for e1, e2 in zip(npos_vec, new_pos_units_list)]
    nvel_vec = BLToCartesian_vel(np.array(pos_vec_vals), np.array(vel_vec_vals), a)
    nvel_vec = [e1 * e2 for e1, e2 in zip(nvel_vec, new_vel_units_list)]
    return (npos_vec, nvel_vec)


def BL2C_8dim(vec, a):
    """
    Function to convert Boyer-Lindquist 8-length numpy array coordinates to Cartesian

    Parameters
    ----------
    vec : ~numpy.array
        Array of shape (n,8) in the form [t,r,theta,phi,vt,vr,vtheta,vphi] in SI units.
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        Array of shape (n,8) in the form [t,x,y,z,vt,vx,vy,vz] in SI units

    """

    newvec = np.copy(vec)

    vectorized_BLToCartesian_pos = np.vectorize(
        BLToCartesian_pos, signature="(n),()->(n)"
    )
    vectorized_BLToCartesian_vel = np.vectorize(
        BLToCartesian_vel, signature="(n),(n),()->(n)"
    )

    newvec[:, 1:4] = vectorized_BLToCartesian_pos(vec[:, 1:4], a)
    newvec[:, 5:8] = vectorized_BLToCartesian_vel(vec[:, 1:4], vec[:, 5:8], a)

    return newvec
