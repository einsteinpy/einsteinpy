import astropy.constants as const
import astropy.units as u
import numpy as np


class Body:
    """
    Class for making a body.

    If the parameters are not passed, they are taken as the parameters of Sun.
    The default body is Sun.

    Parameters
    ----------
    mass - This is the mass of the body.

    pos_vec - This is the position vector of the body. 
              It is a numpy array of 3 elements with their units

    vel_vec - This is the velocity vector of the body. 
              It is a numpy array of 3 elements with their units

    is_attractor - This boolean paramter tells if the body is attractor or not.

    desc - This is the description about the body where the 
           user can describe the object. It defaults to "Empty" string.

    """

    def __init__(
        self,
        mass=const.M_sun.value * u.kg,
        pos_vec=[0 * u.m, 0 * u.rad, 0 * u.rad],
        vel_vec=[0 * u.m / u.s, 0 * u.rad / u.s, 0 * u.rad / u.s],
        is_attractor=False,
        desc="Empty",
    ):

        self.mass = mass
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.is_attractor = is_attractor
        self.desc = desc
