import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45
from einsteinpy.utils import *

_G = constant.G.value
_c = constant.c.value


class Kerr:
    """
    Class for defining Kerr Geometry methods
    """

    @u.quantity.input(time=u.s, M=u.kg, J=(u.kg * (u.m ** 2) / u.s))
    def __init__(self, pos_vec, vel_vec, time, M, J):
        self.M = M
        self.J = J
        self.pos_vec, self.vel_vec = pos_vec, vel_vec
        self.time = time
        # proxy_value
        self.time_vel = 1
        self.initial_vec = np.hstack(
            (self.time.value, self.pos_vec, self.time_vel, self.vel_vec)
        )
        self.schwarzschild_r = schwarzschild_radius(M) * (constant.c ** 2) / constant.G
        self._a = self.J.value / self.M.value

    @classmethod
    def _classmethod_handler(cls, pos_vec, vel_vec, time, M, J):
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
            np.array(pos_vec_vals),
            np.array(vel_vec_vals),
            time.to(u.s),
            M.to(u.kg),
            J.to(u.kg * u.m ** 2 / u.s),
        )

    @classmethod
    def from_BL(cls, pos_vec, vel_vec, time, M, J):
        """
        Constructor

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
        J : ~astropy.units
            Angular Momentum of the body

        """
        cls.input_coord_system = "BL"
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        return cls._classmethod_handler(pos_vec, vel_vec, time, M, J)
