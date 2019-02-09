import astropy.units as u
import numpy as np


def CartesianToSpherical_pos(pos_vec):
    """
    Function to convert cartesian to spherical coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x,y,z coordinates in SI units(m)
    
    Returns
    -------
    a : ~numpy.array
        3-length numpy array with r, theta, phi (m, rad, rad)
    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    tempvar = np.sqrt(pos_vec[0] ** 2 + pos_vec[1] ** 2 + pos_vec[2] ** 2)
    r_vec[0] = tempvar
    r_vec[1] = np.arccos(pos_vec[2] / tempvar)
    r_vec[2] = np.arctan(pos_vec[1] / pos_vec[0])
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
    a : ~numpy.array
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
    a : ~numpy.array
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
    a : ~numpy.array
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
