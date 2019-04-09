import warnings

import astropy.constants as const
import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45
from einsteinpy.utils import *

_G = constant.G.value
_c = constant.c.value


class Body:
    """
    Class for making a General body.

        Parameters
        ----------
        pos_vec : list
            list of 3 components (3 dimensions) along with ~astropy.units
        vel_vec : list
            list of velocities of 3 components (3 dimensions) along with ~astropy.units
        time : ~astropy.units.s
            Time of start
        M : ~astropy.units.kg
            Mass of the body
        a : float
            Spin factor of massive body
        is_attractor : Bool
            To denote is this body is acting as attractor or not
        desc : string
            Description of body which user wants to write        

    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(
        self, pos_vec, vel_vec, M, time=0 * u.s, a=0, is_attractor=False, desc="Empty"
    ):
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.time = time
        self.a = a
        self.M = M
        self.is_attractor = is_attractor
        self.desc = desc

    @classmethod
    def _classmethod_handler(cls, pos_vec, vel_vec, M, time=0 * u.s, a=0):
        cls.units_list = [
            u.s,
            u.m,
            u.rad,
            u.rad,
            u.one,
            u.m / u.s,
            u.rad / u.s,
            u.rad / u.s,
        ]
        pos_vec_vals = [
            pos_vec[i].to(cls.units_list[i + 1]).value for i in range(len(pos_vec))
        ]
        vel_vec_vals = [
            vel_vec[i].to(cls.units_list[i + 5]).value for i in range(len(vel_vec))
        ]
        return cls(
            np.array(pos_vec_vals), np.array(vel_vec_vals), time.to(u.s), M.to(u.kg), a
        )

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_cartesian(cls, pos_vec, vel_vec, M, time=0 * u.s, a=0):
        """
        Constructor
        Initialize the Body object from the Cartesian Coordinates

        Parameters
        ----------
        pos_vec : list
            list of x, y and z components along with ~astropy.units
        vel_vec : list
            list of velocities of x, y, and z components along with ~astropy.units
        time : ~astropy.units.s
            Time of start
        M : ~astropy.units.kg
            Mass of the body
        a : float
            Spin factor of massive body

        """
        cls.input_coord_system = "Cartesian"
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        bl_pos_vec, bl_vel_vec = C2BL_units(pos_vec, vel_vec, a)
        return cls._classmethod_handler(bl_pos_vec, bl_vel_vec, time, M, a)

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_BL(cls, pos_vec, vel_vec, M, time=0 * u.s, a=0):
        """
        Constructor
        Initialize from Boyer-Lindquist Coordinates

        Parameters
        ----------
        pos_vec : list
            list of r, theta & phi components along with ~astropy.units
        vel_vec : list
            list of velocities of r, theta & phi components along with ~astropy.units
        time : ~astropy.units.s
            Time of start
        M : ~astropy.units.kg
            Mass of the body
        a : float
            Spin factor of massive body

        """
        cls.input_coord_system = "Boyer-Lindquist"
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        return cls._classmethod_handler(pos_vec, vel_vec, time, M, a)

    def pos_CatersianToSpherical(self):
        """
        Function to convert the coordinates of position vector from Cartesian to Spherical.
        It used the function CartesianToSperical_pos in coord_transforms module in utils.

        templist1 : List which stores only the values and then converted to numpy array
                    After that it is transformed from Cartesian to Spherical
        templist2 : The values of templist1 are appended to this list and the units are provided.

        """
        templist1 = []
        for i in self.pos_vec:
            templist1.append(i.value)  # get only the values

        templist1 = np.asarray(templist1, dtype=object)  # convert to numpy array
        templist1 = CartesianToSpherical_pos(templist1)  # tranformation

        templist2 = []
        for i in templist1:
            templist2.append(i)

        # Adding units
        templist2[0] = templist2[0] * u.m
        templist2[1] = templist2[1] * u.rad
        templist2[2] = templist2[2] * u.rad
        self.pos_vec = templist2

    def pos_SphericalToCartesian(self):
        """
        Function to convert the coordinates of position vector from Spherical to Cartesian.
        It used the function SphericalToCartesian_pos in coord_transforms module in utils.
        
        templist1 : List which stores only the values and then converted to numpy array
                    After that it is transformed from Spherical to Cartesian.
        templist2 : The values of templist1 are appended to this list and the units are provided.

        """
        templist1 = []
        for i in self.pos_vec:
            templist1.append(i.value)  # get only the values

        templist1 = np.asarray(templist1, dtype=object)  # convert to numpy array
        templist1 = SphericalToCartesian_pos(templist1)  # tranformation

        templist2 = []
        for i in templist1:
            templist2.append(i)

        # Adding units
        templist2[0] = templist2[0] * u.m
        templist2[1] = templist2[1] * u.m
        templist2[2] = templist2[2] * u.m
        self.pos_vec = templist2

    def vel_CartesianToSpherical(self):
        """
        Function to convert the coordinates of velocity vector from Cartesian to Spherical.
        It used the function CartesianToSpherical_vel in coord_transforms module in utils.

        templist1p : List which stores only the values of position vector and then is converted 
                     to numpy array.
        templist1v : List which stores only the values of velocity vector and then is converted 
                     to numpy array. This contains the transformed coordinates from Cartesian to 
                     Spherical.
        templist2v : List which is formed from the numpy array. Units are provided. 

        Note : The position vector must be in the same coordinate system as velocity vector.
        """

        templist1p = []
        for i in self.pos_vec:
            templist1p.append(i.value)  # get only the values
        templist1p = np.asarray(templist1p, dtype=object)  # convert to numpy array

        templist1v = []
        for i in self.vel_vec:
            templist1v.append(i.value)  # get only the values
        templist1v = np.asarray(templist1v, dtype=object)  # convert to numpy array

        templist1v = CartesianToSpherical_vel(templist1p, templist1v)  # tranformation

        templist2v = []
        for i in templist1v:
            templist2v.append(i)

        # Adding units
        templist2v[0] = templist2v[0] * u.m / u.s
        templist2v[1] = templist2v[1] * u.rad / u.s
        templist2v[2] = templist2v[2] * u.rad / u.s
        self.vel_vec = templist2v

    def vel_SphericalToCartesian(self):
        """
        Function to convert the coordinates of velocity vector from Spherical to Cartesian.
        It used the function SphericalToCartesian_vel in coord_transforms module in utils.
        
        templist1p : List which stores only the values of position vector and then is converted 
                     to numpy array.
        templist1v : List which stores only the values of velocity vector and then is converted 
                     to numpy array. This contains the transformed coordinates from Spherical to 
                     Cartesian.
        templist2v : List which is formed from the numpy array. Units are provided.

        Note : The position vector must be in the same coordinate system as velocity vector. 

        """
        templist1p = []
        for i in self.pos_vec:
            templist1p.append(i.value)  # get only the values
        templist1p = np.asarray(templist1p, dtype=object)  # convert to numpy array

        templist1v = []
        for i in self.vel_vec:
            templist1v.append(i.value)  # get only the values
        templist1v = np.asarray(templist1v, dtype=object)  # convert to numpy array

        templist1v = SphericalToCartesian_vel(templist1p, templist1v)  # tranformation

        templist2v = []
        for i in templist1v:
            templist2v.append(i)

        # Adding units
        templist2v[0] = templist2v[0] * u.m / u.s
        templist2v[1] = templist2v[1] * u.m / u.s
        templist2v[2] = templist2v[2] * u.m / u.s
        self.vel_vec = templist2v
