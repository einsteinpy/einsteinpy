import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.coordinates import BoyerLindquistDifferential, CartesianDifferential
from einsteinpy.integrators import RK45
from einsteinpy.utils import kerrnewman_utils, schwarzschild_radius

_G = constant.G.value
_c = constant.c.value
_Cc = constant.coulombs_const.value
# q = 11604461683.91822052001953125*q This comment is not to be removed. This is a hack.


class KerrNewman:
    """
    Class for defining Kerr-Newman Goemetry Methods
    """

    @u.quantity_input(time=u.s, M=u.kg, Q=u.C)
    def __init__(self, bl_coords, q, M, Q, time):
        self.input_coords = bl_coords
        self.M = M
        self.a = self.input_coords.a.to(u.m)
        self.Q = Q
        self.q = q

        pos_vec, vel_vec = (
            self.input_coords.si_values()[:3],
            self.input_coords.si_values()[3:],
        )
        self.time = time
        time_vel = kerrnewman_utils.kerrnewman_time_velocity(
            pos_vec, vel_vec, self.M, self.a.value, Q
        )
        self.initial_vec = np.hstack(
            (self.time.value, pos_vec, time_vel.value, vel_vec)
        )
        self.scr = schwarzschild_radius(M)

    @classmethod
    @u.quantity_input(q=u.C / u.kg, time=u.s, M=u.kg, Q=u.C, a=u.m)
    def from_coords(cls, coords, M, q, Q, time=0 * u.s, a=0 * u.m):
        """
        Constructor.

        Parameters
        ----------
        coords : ~einsteinpy.coordinates.velocity.BoyerLindquistDifferential
            Initial positions and velocities of particle in BL Coordinates, and spin factor of massive body.
        q : ~astropy.units.quantity.Quantity
            Charge per unit mass of test particle
        M : ~astropy.units.quantity.Quantity
            Mass of the massive body
        Q : ~astropy.units.quantity.Quantity
            Charge on the massive body
        a : ~astropy.units.quantity.Quantity
            Spin factor of the massive body(Angular Momentum per unit mass per speed of light)
        time : ~astropy.units.quantity.Quantity
            Time of start, defaults to 0 seconds.

        """
        if coords.system == "Cartesian":
            bl_coords = coords.bl_differential(a)
            return cls(bl_coords, q, M, Q, time)
        if coords.system == "Spherical":
            bl_coords = coords.bl_differential(a)
            return cls(bl_coords, q, M, Q, time)
        return cls(coords, q, M, Q, time)

    def f_vec(self, ld, vec):
        chl = kerrnewman_utils.christoffels(
            vec[1], vec[2], self.M.value, self.a.value, self.Q.value
        )
        maxwell = kerrnewman_utils.maxwell_tensor_contravariant(
            vec[1], vec[2], self.a.value, self.Q.value, self.M.value
        )
        metric = kerrnewman_utils.metric(
            vec[1], vec[2], self.M.value, self.a.value, self.Q.value
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
        vals[4:] -= self.q.value * np.dot(
            vec[4:].reshape((4,)), np.matmul(metric, maxwell)
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
            Starting lambda(proper time), defaults to 0.0, (lambda ~= t)
        end_lamdba : float
            Lambda(proper time) where iteartions will stop, defaults to 100000
        stop_on_singularity : bool
            Whether to stop further computation on reaching singularity, defaults to True
        OdeMethodKwargs : dict
            Kwargs to be supplied to the ODESolver, defaults to {'stepsize': 1e-3}
            Dictionary with key 'stepsize' along with an float value is expected.
        return_cartesian : bool
            True if coordinates and velocities are required in cartesian coordinates(SI units), defaults to False

        Returns
        -------
        ~numpy.ndarray
            N-element array containing proper time.
        ~numpy.ndarray
            (n,8) shape array of [t, x1, x2, x3, velocity_of_time, v1, v2, v3] for each proper time(lambda).

        """
        vecs = list()
        lambdas = list()
        crossed_event_horizon = False
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs
        )

        _scr = self.scr.value * 1.001

        while ODE.t < end_lambda:
            vecs.append(ODE.y)
            lambdas.append(ODE.t)
            ODE.step()
            if (not crossed_event_horizon) and (ODE.y[1] <= _scr):
                warnings.warn("particle reached Schwarzchild Radius. ", RuntimeWarning)
                if stop_on_singularity:
                    break
                else:
                    crossed_event_horizon = True

        vecs, lambdas = np.array(vecs), np.array(lambdas)

        if not return_cartesian:
            return lambdas, vecs
        else:
            cart_vecs = list()
            for v in vecs:
                si_vals = (
                    BoyerLindquistDifferential(
                        v[1] * u.m,
                        v[2] * u.rad,
                        v[3] * u.rad,
                        v[5] * u.m / u.s,
                        v[6] * u.rad / u.s,
                        v[7] * u.rad / u.s,
                        self.a,
                    )
                    .cartesian_differential()
                    .si_values()
                )
                cart_vecs.append(np.hstack((v[0], si_vals[:3], v[4], si_vals[3:])))
            return lambdas, np.array(cart_vecs)

    def calculate_trajectory_iterator(
        self,
        start_lambda=0.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation.
        Yields an iterator.

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
        float
            proper time
        ~numpy.ndarray
            array of [t, x1, x2, x3, velocity_of_time, v1, v2, v3] for each proper time(lambda).

        """
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=1e300,
            **OdeMethodKwargs
        )
        crossed_event_horizon = False
        _scr = self.scr.value * 1.001

        while True:
            if not return_cartesian:
                yield ODE.t, ODE.y
            else:
                v = ODE.y
                si_vals = (
                    BoyerLindquistDifferential(
                        v[1] * u.m,
                        v[2] * u.rad,
                        v[3] * u.rad,
                        v[5] * u.m / u.s,
                        v[6] * u.rad / u.s,
                        v[7] * u.rad / u.s,
                        self.a,
                    )
                    .cartesian_differential()
                    .si_values()
                )
                yield ODE.t, np.hstack((v[0], si_vals[:3], v[4], si_vals[3:]))
            ODE.step()
            if (not crossed_event_horizon) and (ODE.y[1] <= _scr):
                warnings.warn("particle reached Schwarzschild Radius. ", RuntimeWarning)
                if stop_on_singularity:
                    break
                else:
                    crossed_event_horizon = True
