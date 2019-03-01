import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45, RK4naive
from einsteinpy.utils import *
from einsteinpy.utils import schwarzschild_radius as scr
from einsteinpy.utils import time_velocity

_G = constant.G.value
_c = constant.c.value


class Schwarzschild:
    """
    Class for defining a Schwarzschild Geometry methods
    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(self, pos_vec, vel_vec, time, M):
        self.M = M
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.time = time
        self.time_vel = time_velocity(pos_vec, vel_vec, M)
        self.initial_vec = np.hstack(
            (self.time.value * _c, self.pos_vec, self.time_vel.value, self.vel_vec / _c)
        )
        self.schwarzschild_r = scr(M)

    @classmethod
    def _classmethod_handler(cls, pos_vec, vel_vec, time, M):
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
            np.array(pos_vec_vals), np.array(vel_vec_vals), time.to(u.s), M.to(u.kg)
        )

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_spherical(cls, pos_vec, vel_vec, time, M):
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

        """
        cls.input_coord_system = "Spherical"
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        return cls._classmethod_handler(pos_vec, vel_vec, time, M)

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg)
    def from_cartesian(cls, pos_vec, vel_vec, time, M):
        """
        Constructor

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

        """
        cls.input_coord_system = "Cartesian"
        cls.input_units_list = (
            [time.unit]
            + [pos_vec[i].unit for i in range(len(pos_vec))]
            + [u.one]
            + [vel_vec[i].unit for i in range(len(vel_vec))]
        )
        sp_pos_vec, sp_vel_vec = C2S_units(pos_vec, vel_vec)
        return cls._classmethod_handler(sp_pos_vec, sp_vel_vec, time, M)

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
        def f0_3():
            return vec[i + 4]

        def f4():
            term1 = self.christ_sym0_01(vec) * vec[4] * vec[5]
            return -2 * term1

        def f5():
            term1 = self.christ_sym1_00(vec) * vec[4] * vec[4]
            term2 = self.christ_sym1_11(vec) * vec[5] * vec[5]
            term3 = self.christ_sym1_22(vec) * vec[6] * vec[6]
            term4 = self.christ_sym1_33(vec) * vec[7] * vec[7]
            return -1 * (term1 + term2 + term3 + term4)

        def f6():
            term1 = self.christ_sym2_21(vec) * vec[6] * vec[5]
            term2 = self.christ_sym2_33(vec) * vec[7] * vec[7]
            return -1 * (2 * term1 + term2)

        def f7():
            term1 = self.christ_sym3_31(vec) * vec[7] * vec[5]
            term2 = self.christ_sym3_32(vec) * vec[7] * vec[6]
            return -1 * (2 * term1 + 2 * term2)

        f_dict = {0: f0_3, 1: f0_3, 2: f0_3, 3: f0_3, 4: f4, 5: f5, 6: f6, 7: f7}
        return f_dict[i]()

    def f_vec(self, ld, vec):
        f_vec_vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
        for i in range(len(vec)):
            f_vec_vals[i] = self.f(i, vec)
        return _c * f_vec_vals

    def calculate_trajectory(
        self,
        start_lambda=0.0,
        end_lambda=10.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation

        Parameters
        ----------
        start_lambda : float
            Starting lambda, defaults to 0.0, (lambda ~= t)
        end_lamdba : float
            Lambda where iteartions will stop, defaults to 100000
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity, defaults to True
        OdeMethodKwargs : dict
            Kwargs to be supplied to the ODESolver, defaults to {'stepsize': 1e-3}
            Dictionary with key 'stepsize' along with an float value is expected.
        return_cartesian : bool
            True if coordinates and velocities are required in cartesian coordinates(SI units), defaults to False

        Returns
        -------
        tuple
            (~numpy.array of lambda, (n,8) shape ~numpy.array of [t, pos1, pos2, pos3, velocity_of_time, vel1, vel2, vel3])

        """
        vec_list = list()
        lambda_list = list()
        singularity_reached = False
        scaling_factors = np.array([1 / _c, 1.0, 1.0, 1.0, 1.0, _c, _c, _c])
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs
        )
        _scr = self.schwarzschild_r.value * 1.001
        while ODE.t < end_lambda:
            vec_list.append(ODE.y)
            lambda_list.append(ODE.t)
            ODE.step()
            if (not singularity_reached) and (ODE.y[1] <= _scr):
                warnings.warn(
                    "r component of position vector reached Schwarzchild Radius. ",
                    RuntimeWarning,
                )
                if stop_on_singularity:
                    break
                else:
                    singularity_reached = True

        def _not_cartesian():
            return (np.array(lambda_list), np.array(vec_list) * scaling_factors)

        def _cartesian():
            self.units_list = [
                u.s,
                u.m,
                u.m,
                u.m,
                u.one,
                u.m / u.s,
                u.m / u.s,
                u.m / u.s,
            ]
            return (
                np.array(lambda_list),
                S2C_8dim(np.array(vec_list) * scaling_factors),
            )

        choice_dict = {0: _not_cartesian, 1: _cartesian}
        return choice_dict[int(return_cartesian)]()

    def calculate_trajectory_iterator(
        self,
        start_lambda=0.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation
        Yields an iterator

        Parameters
        ----------
        start_lambda : float
            Starting lambda, defaults to 0.0, (lambda ~= t)
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity, defaults to True
        OdeMethodKwargs : dict
            Kwargs to be supplied to the ODESolver, defaults to {'stepsize': 1e-3}
            Dictionary with key 'stepsize' along with an float value is expected.
        return_cartesian : bool
            True if coordinates and velocities are required in cartesian coordinates(SI units), defaults to Falsed

        Yields
        ------
        tuple
            (lambda, (8,) shape ~numpy.array of [t, pos1, pos2, pos3, velocity_of_time, vel1, vel2, vel3])

        """
        singularity_reached = False
        scaling_factors = np.array([1 / _c, 1.0, 1.0, 1.0, 1.0, _c, _c, _c])
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=1e300,
            **OdeMethodKwargs
        )
        _scr = self.schwarzschild_r.value * 1.001

        def yielder_func():
            nonlocal singularity_reached
            while True:
                if not return_cartesian:
                    yield (ODE.t, np.multiply(ODE.y, scaling_factors))
                else:
                    temp = np.copy(ODE.y)
                    temp[0] *= scaling_factors[0]
                    temp[1:4] = SphericalToCartesian_pos(
                        np.multiply(ODE.y[1:4], scaling_factors[1:4])
                    )
                    temp[5:8] = SphericalToCartesian_vel(
                        np.multiply(ODE.y[1:4], scaling_factors[1:4]),
                        np.multiply(ODE.y[5:8], scaling_factors[5:8]),
                    )
                    yield (ODE.t, temp)
                ODE.step()
                if (not singularity_reached) and (ODE.y[1] <= _scr):
                    warnings.warn(
                        "r component of position vector reached Schwarzchild Radius. ",
                        RuntimeWarning,
                    )
                    if stop_on_singularity:
                        break
                    else:
                        singularity_reached = True

        if return_cartesian:
            self.units_list = [
                u.s,
                u.m,
                u.m,
                u.m,
                u.one,
                u.m / u.s,
                u.m / u.s,
                u.m / u.s,
            ]
        return yielder_func()
