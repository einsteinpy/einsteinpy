import warnings

import astropy.units as u
import numpy as np

warnings.warn(
    "This module will be deprecated in EinsteinPy v0.3.0. Please switch to `einsteinpy.coordinates`.",
    PendingDeprecationWarning,
)


def CartesianToSpherical_pos(pos_vec):
    """
    Function to convert cartesian to spherical coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x,y,z coordinates in SI units(m)

    Returns
    -------
    ~numpy.array
        3-length numpy array with r, theta, phi (m, rad, rad)
    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    tempvar = np.sqrt(pos_vec[0] ** 2 + pos_vec[1] ** 2 + pos_vec[2] ** 2)
    r_vec[0] = tempvar
    r_vec[1] = np.arccos(pos_vec[2] / tempvar)
    r_vec[2] = np.arctan2(pos_vec[1], pos_vec[0])
    return r_vec


def CartesianToSpherical_vel(pos_vec, vel_vec):
    """
    Function to convert velocities in cartesian coordinates to velocities in spherical coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x,y,z coordinates in SI units(m)
    vel_vec : ~numpy.array
        3-length numpy array having vx, vy, vz in SI units(m/s)

    Returns
    -------
    ~numpy.array
        3-length numpy array with V_r, V_theta, V_phi (m/s, rad/s, rad/s)

    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    tempvar1 = pos_vec[0] ** 2 + pos_vec[1] ** 2
    tempvar2 = tempvar1 + pos_vec[2] ** 2
    v_vec[0] = (
        pos_vec[0] * vel_vec[0] + pos_vec[1] * vel_vec[1] + pos_vec[2] * vel_vec[2]
    ) / np.sqrt(tempvar2)
    v_vec[2] = -1 * (vel_vec[0] * pos_vec[1] - pos_vec[0] * vel_vec[1]) / tempvar1
    v_vec[1] = (
        pos_vec[2] * (pos_vec[0] * vel_vec[0] + pos_vec[1] * vel_vec[1])
        - tempvar1 * vel_vec[2]
    ) / (tempvar2 * np.sqrt(tempvar1))
    return v_vec


def SphericalToCartesian_pos(pos_vec):
    """
    Function to convert spherical coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)

    Returns
    -------
    ~numpy.array
        3-length numpy array with x, y, z in m.

    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    r_vec[0] = pos_vec[0] * np.cos(pos_vec[2]) * np.sin(pos_vec[1])
    r_vec[1] = pos_vec[0] * np.sin(pos_vec[2]) * np.sin(pos_vec[1])
    r_vec[2] = pos_vec[0] * np.cos(pos_vec[1])
    return r_vec


def SphericalToCartesian_vel(pos_vec, vel_vec):
    """
    Function to convert velocities in spherical coordinates to velocities in cartesian coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    vel_vec : ~numpy.array
        3-length numpy array with V_r, V_theta, V_phi (m/s, rad/s, rad/s)

    Returns
    -------
    ~numpy.array
        3-length numpy array having vx, vy, vz in SI units(m/s)

    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    v_vec[0] = (
        np.sin(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[0]
        - pos_vec[0] * np.sin(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[2]
        + pos_vec[0] * np.cos(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[1]
    )
    v_vec[1] = (
        np.sin(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[0]
        + pos_vec[0] * np.cos(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[1]
        + pos_vec[0] * np.sin(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[2]
    )
    v_vec[2] = (
        np.cos(pos_vec[1]) * vel_vec[0] - pos_vec[0] * np.sin(pos_vec[1]) * vel_vec[1]
    )
    return v_vec


def C2S_units(pos_vec, vel_vec):
    """
    Function to convert Cartesian to Spherical Coordinates along with handling units

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
    npos_vec = CartesianToSpherical_pos(np.array(pos_vec_vals)).tolist()
    npos_vec = [e1 * e2 for e1, e2 in zip(npos_vec, new_pos_units_list)]
    nvel_vec = CartesianToSpherical_vel(np.array(pos_vec_vals), np.array(vel_vec_vals))
    nvel_vec = [e1 * e2 for e1, e2 in zip(nvel_vec, new_vel_units_list)]
    return (npos_vec, nvel_vec)


def S2C_units(pos_vec, vel_vec):
    """
    Function to convert Spherical to Cartesian Coordinates along with handling units

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
    npos_vec = SphericalToCartesian_pos(np.array(pos_vec_vals)).tolist()
    npos_vec = [e1 * e2 for e1, e2 in zip(npos_vec, new_pos_units_list)]
    nvel_vec = SphericalToCartesian_vel(np.array(pos_vec_vals), np.array(vel_vec_vals))
    nvel_vec = [e1 * e2 for e1, e2 in zip(nvel_vec, new_vel_units_list)]
    return (npos_vec, nvel_vec)


def S2C_8dim(vec):
    """
    Function to convert spherical 8-length numpy array coordinates to cartesian

    Parameters
    ----------
    vec : ~numpy.array
        Array of shape (n,8) in the form [t,r,theta,phi,vt,vr,vtheta,vphi] in SI units.

    Returns
    -------
    ~numpy.array
        Array of shape (n,8) in the form [t,x,y,z,vt,vx,vy,vz] in SI units

    """

    newvec = np.copy(vec)

    vectorized_SphericalToCartesian_pos = np.vectorize(
        SphericalToCartesian_pos, signature="(n)->(n)"
    )
    vectorized_SphericalToCartesian_vel = np.vectorize(
        SphericalToCartesian_vel, signature="(n),(n)->(n)"
    )

    newvec[:, 1:4] = vectorized_SphericalToCartesian_pos(vec[:, 1:4])
    newvec[:, 5:8] = vectorized_SphericalToCartesian_vel(vec[:, 1:4], vec[:, 5:8])

    return newvec
