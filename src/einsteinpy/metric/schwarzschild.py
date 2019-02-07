import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45, RK4naive
from einsteinpy.utils import schwarzschild_radius, time_velocity

_G = constant.G.value
_c = constant.c.value


class Schwarzschild:
    """
    Class for defining a Schwarzschild Metric
    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(self, pos_vec, vel_vec, time, M):
        self.M = M
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.time = time
        self.time_vel = time_velocity(pos_vec, vel_vec, M)
        self.initial_vec = np.hstack(
            (time.value * _c, self.pos_vec, self.time_vel.value, self.vel_vec / _c)
        )
        self.schwarzschild_r = schwarzschild_radius(M)

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_spherical(cls, pos_vec, vel_vec, time, M):
        """
        Constructor

        Parameters
        ----------
        pos_vector : list
            list of r, theta & phi components along with astropy units
        vel_vector : list
            list of velocities of r, theta & phi components along with astropy units
        time : ~astropy.units.s
            Time of start
        M : ~astropy.units.kg
            Mass of the body

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
            np.array(pos_vec_vals), np.array(vel_vec_vals), time.to(u.s), M.to(u.kg)
        )

    def christ_sym1_00(self, vec):
        num1 = (-2 * _G * self.M.value) + ((_c ** 2) * vec[1])
        num2 = _G * self.M.value
        deno1 = _c ** 4
        deno2 = vec[1] ** 3
        return (num1 / deno1) * (num2 / deno2)

    def christ_sym1_11(self, vec):
        num = _G * self.M.value
        deno1 = 2 * _G * self.M.value * vec[1]
        deno2 = (_c ** 2) * (vec[1] ** 2)
        deno = deno1 - deno2
        return num / deno

    def christ_sym1_22(self, vec):
        return self.schwarzschild_r.value - vec[1]

    def christ_sym1_33(self, vec):
        return (
            (-1 + (self.schwarzschild_r.value / vec[1])) * vec[1] * np.sin(vec[2]) ** 2
        )

    def christ_sym0_01(self, vec):
        num = _G * self.M.value
        deno1 = -2 * _G * self.M.value + (_c ** 2) * vec[1]
        deno2 = vec[1]
        return (num / deno1) * (1 / deno2)

    def christ_sym2_21(self, vec):
        return 1 / vec[1]

    def christ_sym2_33(self, vec):
        return -1 * np.cos(vec[2]) * np.sin(vec[2])

    def christ_sym3_31(self, vec):
        return 1 / vec[1]

    def christ_sym3_32(self, vec):
        return np.cos(vec[2]) / np.sin(vec[2])

    def f(self, i, vec):
        if i in range(4):
            value = vec[i + 4]
        elif i == 4:
            term1 = self.christ_sym0_01(vec) * vec[4] * vec[5]
            value = -2 * term1
        elif i == 5:
            term1 = self.christ_sym1_00(vec) * vec[4] * vec[4]
            term2 = self.christ_sym1_11(vec) * vec[5] * vec[5]
            term3 = self.christ_sym1_22(vec) * vec[6] * vec[6]
            term4 = self.christ_sym1_33(vec) * vec[7] * vec[7]
            value = -1 * (term1 + term2 + term3 + term4)
        elif i == 6:
            term1 = self.christ_sym2_21(vec) * vec[6] * vec[5]
            term2 = self.christ_sym2_33(vec) * vec[7] * vec[7]
            value = -1 * (2 * term1 + term2)
        elif i == 7:
            term1 = self.christ_sym3_31(vec) * vec[7] * vec[5]
            term2 = self.christ_sym3_32(vec) * vec[7] * vec[6]
            value = -1 * (2 * term1 + term2)

        return value

    def f_vec(self, ld, vec):
        f_vec_vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
        for i in range(len(vec)):
            f_vec_vals[i] = self.f(i, vec)
        return _c * f_vec_vals

    def calculate_trajectory(
        self,
        start_lambda=0.0,
        end_lambda=1e7,
        stop_on_singularity=True,
        OdeMethodKwargs={},
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
        singularity_reached = False
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
            if (not singularity_reached) and (ODE.y[1] <= self.schwarzschild_r.value):
                warnings.warn(
                    "r component of position vector reached Schwarzchild Radius. ",
                    RuntimeWarning,
                )
                if stop_on_singularity:
                    break
                else:
                    singularity_reached = True
        scaling_factors = np.array([1 / _c, 1.0, 1.0, 1.0, 1.0, _c, _c, _c])
        return (np.array(lambda_list), np.array(vec_list) * scaling_factors)
