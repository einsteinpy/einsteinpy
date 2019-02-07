import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45, RK4naive
from einsteinpy.utils import scalar_factor as sf, scalar_factor_derivative as sfd
from einsteinpy.utils import time_velocity

_G = constant.G.value
_c = constant.c.value


class RobertsonWalker:
    """
    Class for defining a Robertson-Walker Geometry methods
    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(self, pos_vec, vel_vec, time, M, era, tuning_param, k):
        self.M = M
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.time = time
        self.time_vel = time_velocity(pos_vec, vel_vec, M)
        self.initial_vec = np.hstack(
            (time.value, self.pos_vec, self.time_vel.value, self.vel_vec)
        )
        self.era = era
        self.tuning_param = tuning_param
        self.k = k

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_spherical(
        cls, pos_vec, vel_vec, time, M, era="md", tuning_param=1.0, k=0.0
    ):
        """
        Constructor

        Parameters
        ----------
        pos_vec : list
            list of r, theta & phi components along with astropy units
        vel_vec : list
            list of velocities of r, theta & phi components along with astropy units
        time : ~astropy.units.s
            Time of start
        M : ~astropy.units.kg
            Mass of the body
        era : string
            Can be chosen from 'md' (Matter Dominant),
            'rd' (Radiation Dominant) and 'ded' (Dark Energy Dominant)
        tuning_param : float, optional
            Unit scaling factor, defaults to 1
        k : float, optional
            Curvature of universe, 0 for flat, 1 for postive curvature, -1 for negetive curvature

        """
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
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        return cls(
            np.array(pos_vec_vals),
            np.array(vel_vec_vals),
            time.to(u.s),
            M.to(u.kg),
            era,
            tuning_param,
            k,
        )

    def christ_sym0_11(self, vec):
        num = sf(vec[0] * u.s, self.era, self.tuning_param) * sfd(
            vec[0] * u.s, self.era, self.tuning_param
        )
        deno = 1 - (self.k * (vec[1] ** 2))
        return num / deno

    def christ_sym0_22(self, vec):
        pass
