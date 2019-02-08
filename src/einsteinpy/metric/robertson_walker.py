import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45, RK4naive
from einsteinpy.utils import scalar_factor as sf
from einsteinpy.utils import scalar_factor_derivative as sfd
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
        val = (
            sf(vec[0] * u.s, self.era, self.tuning_param)
            * sfd(vec[0] * u.s, self.era, self.tuning_param)
            * (vec[1] ** 2)
        )
        return val

    def christ_sym0_33(self, vec):
        val = (
            sf(vec[0] * u.s, self.era, self.tuning_param)
            * sfd(vec[0] * u.s, self.era, self.tuning_param)
            * (vec[1] ** 2)
            * (np.sin(vec[2]) ** 2)
        )
        return val

    def christ_sym1_01(self, vec):
        val = sfd(vec[0] * u.s, self.era, self.tuning_param) / sf(
            vec[0] * u.s, self.era, self.tuning_param
        )
        return val

    def christ_sym1_10(self, vec):
        return self.christ_sym1_01(vec)

    def christ_sym2_02(self, vec):
        return self.christ_sym1_01(vec)

    def christ_sym2_20(self, vec):
        return self.christ_sym1_01(vec)

    def christ_sym3_03(self, vec):
        return self.christ_sym1_01(vec)

    def christ_sym3_30(self, vec):
        return self.christ_sym1_01(vec)

    def christ_sym2_12(self, vec):
        val = 1.0 / vec[1]
        return val

    def christ_sym2_21(self, vec):
        return self.christ_sym2_12(vec)

    def christ_sym3_13(self, vec):
        return self.christ_sym2_12(vec)

    def christ_sym3_31(self, vec):
        return self.christ_sym2_12(vec)

    def christ_sym2_33(self, vec):
        val = -1.0 * np.sin(vec[2]) * np.cos(vec[2])
        return val

    def christ_sym3_23(self, vec):
        val = np.cos(vec[2]) / np.sin(vec[2])
        return val

    def christ_sym3_32(self, vec):
        return self.christ_sym3_23(vec)

    def f(self, i, vec):
        def f0_3():
            return vec[i + 4]

        def f4():
            term1 = self.christ_sym0_11(vec) * vec[5] * vec[5]
            term2 = self.christ_sym0_22(vec) * vec[6] * vec[6]
            term3 = self.christ_sym0_33(vec) * vec[7] * vec[7]
            return -1 * (term1 + term2 + term3)

        def f5():
            term1 = self.christ_sym1_01(vec) * vec[4] * vec[5]
            return -2 * term1

        def f6():
            term1 = self.christ_sym2_12(vec) * vec[5] * vec[6]
            term2 = self.christ_sym2_33(vec) * vec[7] * vec[7]
            return -1 * (2 * term1 + term2)

        def f7():
            term1 = self.christ_sym3_31(vec) * vec[7] * vec[5]
            term2 = self.christ_sym3_23(vec) * vec[6] * vec[7]
            return -2 * (term1 + term2)

        f_dict = {0: f0_3, 1: f0_3, 2: f0_3, 3: f0_3, 4: f4, 5: f5, 6: f6, 7: f7}
        return f_dict[i]()

    def f_vec(self, ld, vec):
        f_vec_vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
        for i in range(len(vec)):
            f_vec_vals[i] = self.f(i, vec)
        return f_vec_vals

    def calculate_trajectory(
        self, start_lambda=0.0, end_lambda=10.0, OdeMethodKwargs={"stepsize": 1e-3}
    ):
        """
        Calculate trajectory in manifold according to geodesic equation

        Parameters
        ----------
        start_lambda : float
            Starting lambda, defaults to 0.0, ( c ~= t)
        end_lamdba : float
            Lambda where iteartions will stop, defaults to 100000
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity, defaults to True

        Returns
        -------
        tuple of 2 numpy.array

        """
        vec_list = list()
        lambda_list = list()
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs
        )
        while ODE.t < end_lambda:
            vec_list.append(ODE.y)
            lambda_list.append(ODE.t)
            ODE.step()
        scaling_factors = 1.0
        # scaling_factors = np.array([1 / _c, 1.0, 1.0, 1.0, 1.0, _c, _c, _c])
        return (np.array(lambda_list), np.array(vec_list) * scaling_factors)
