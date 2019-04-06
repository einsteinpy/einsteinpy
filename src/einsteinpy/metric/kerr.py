import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.integrators import RK45
from einsteinpy.utils import *
from einsteinpy.utils import kerr_utils
from einsteinpy.utils import schwarzschild_radius as scr

_G = constant.G.value
_c = constant.c.value


class Kerr:
    """
    Class for defining Kerr Goemetry Methdos
    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(self, pos_vec, vel_vec, time, M, a):
        self.M = M
        self.a = a
        self.pos_vec = pos_vec
        self.vel_vec = vel_vec
        self.time = time
        self.time_vel = kerr_utils.kerr_time_velocity(
            self.pos_vec, self.vel_vec, self.M, self.a
        )
        self.initial_vec = np.hstack(
            (self.time.value, self.pos_vec, self.time_vel.value, self.vel_vec)
        )
        self.schwarzschild_r = scr(M)
        self.kerr_r = kerr_utils.event_horizon(self.schwarzschild_r, self.a, theta=np.pi / 2)

    @classmethod
    def _classmethod_handler(cls, pos_vec, vel_vec, time, M, a):
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
    def from_BL(cls, pos_vec, vel_vec, time, M, a):
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

    @classmethod
    def from_cartesian(cls, pos_vec, vel_vec, time, M, a):
        """
        Constructor
        Initialize from Cartesian Coordinates

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

    def f_vec(self, ld, vec):
        chl = kerr_utils.christoffels(
            _c, vec[1], vec[2], self.schwarzschild_r.value, self.a
        )
        vals = np.zeros(shape=(8,), dtype=float)
        for i in range(4):
            vals[i] = vec[i + 4]
        vals[4] = -2.0 * (
            chl[0][0][1] * vec[4] * vec[5]
            + chl[0][0][2] * vec[4] * vec[6]
            + chl[0][1][3] * vec[5] * vec[7]
            + chl[0][2][3] * vec[6] * vec[7]
        )
        vals[5] = -1.0 * (
            chl[1][0][0] * vec[4] * vec[4]
            + 2 * chl[1][0][3] * vec[4] * vec[7]
            + chl[1][1][1] * vec[5] * vec[5]
            + 2 * chl[1][1][2] * vec[5] * vec[6]
            + chl[1][2][2] * vec[6] * vec[6]
            + chl[1][3][3] * vec[7] * vec[7]
        )
        vals[6] = -1.0 * (
            chl[2][0][0] * vec[4] * vec[4]
            + 2 * chl[2][0][3] * vec[4] * vec[7]
            + chl[2][1][1] * vec[5] * vec[5]
            + 2 * chl[2][1][2] * vec[5] * vec[6]
            + chl[2][2][2] * vec[6] * vec[6]
            + chl[2][3][3] * vec[7] * vec[7]
        )
        vals[7] = -2.0 * (
            chl[3][0][1] * vec[4] * vec[5]
            + chl[3][0][2] * vec[4] * vec[6]
            + chl[3][1][3] * vec[5] * vec[7]
            + chl[3][2][3] * vec[6] * vec[7]
        )
        return vals

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
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs
        )
        _scr = self.schwarzschild_r.value * 1.001
        _kerr_r = kerr_r * 1.001
        while ODE.t < end_lambda:
            vec_list.append(ODE.y)
            lambda_list.append(ODE.t)
            ODE.step()
            if (not singularity_reached) and (ODE.y[1] <= _kerr_r):
                warnings.warn(
                    "r component of position vector reached Kerr Radius. ",
                    RuntimeWarning,
                )
                if stop_on_singularity:
                    break
                else:
                    singularity_reached = True

        def _not_cartesian():
            return (np.array(lambda_list), np.array(vec_list))

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
            return (np.array(lambda_list), BL2C_8dim(np.array(vec_list), self.a))

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
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=1e300,
            **OdeMethodKwargs
        )
        _scr = self.schwarzschild_r.value * 1.001
        _kerr_r = kerr*1.001

        def yielder_func():
            nonlocal singularity_reached
            while True:
                if not return_cartesian:
                    yield (ODE.t, ODE.y)
                else:
                    temp = np.copy(ODE.y)
                    temp[1:4] = BLToCartesian_pos(ODE.y[1:4], self.a)
                    temp[5:8] = BLToCartesian_vel(ODE.y[1:4], ODE.y[5:8], self.a)
                    yield (ODE.t, temp)
                ODE.step()
                if (not singularity_reached) and (ODE.y[1] <= _kerr_r):
                    warnings.warn(
                        "r component of position vector reached Kerr Radius. ",
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
