import astropy.units as u
import numpy as np

from einsteinpy import constant


@u.quantity_input(mass=u.kg)
def schwarzschild_radius(mass):
    """
    Schwarzschild radius is the radius defining the event horizon of a
    Schwarzschild black hole. It is characteristic radius associated with every
    quantity of mass.

    Parameters
    ----------
    mass : ~astropy.units.kg

    Returns
    -------
    r : ~astropy.units.m
        Schwarzschild radius for a given mass

    """
    M = mass.to(u.kg)
    num = 2 * constant.G * M
    deno = constant.c ** 2
    return num / deno


@u.quantity_input(mass=u.kg)
def time_velocity(pos_vec, vel_vec, mass):
    """
    Velocity of time calculated from einstein's equation.
    See http://www.physics.usu.edu/Wheeler/GenRel/Lectures/GRNotesDecSchwarzschildGeodesicsPost.pdf

    Parameters
    ----------
    pos_vector : ~numpy.array
        Vector with r, theta, phi components in SI units
    vel_vector : ~numpy.array
        Vector with velocities of r, theta, phi components in SI units
    mass : ~astropy.units.kg
        Mass of the body

    Returns
    -------
    ~astropy.units.one
        Velocity of time

    """
    # this function considers SI units only
    a = schwarzschild_radius(mass).value
    c = constant.c.value
    num1 = (1 / ((c ** 2) * (1 + (a / pos_vec[0])))) * (vel_vec[0] ** 2)
    num2 = ((pos_vec[0] ** 2) / (c ** 2)) * (vel_vec[1] ** 2)
    num3 = (
        ((pos_vec[0] ** 2) / (c ** 2)) * (np.sin(pos_vec[1]) ** 2) * (vel_vec[2] ** 2)
    )
    deno = 1 + (a / pos_vec[0])
    time_vel_squared = (1 + num1 + num2 + num3) / deno
    time_vel = np.sqrt(time_vel_squared)
    return time_vel * u.one
