import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.coordinates import CartesianDifferential, SphericalDifferential
from einsteinpy.integrators import RK45
from einsteinpy.utils import schwarzschild_radius, schwarzschild_utils

_G = constant.G.value
_c = constant.c.value


class Schwarzschild:
    """
    Class for defining a Schwarzschild Geometry methods
    """

    @u.quantity_input(time=u.s, M=u.kg)
    def __init__(self, sph_coords, M, time):
        self.M = M
        self.input_coords = sph_coords
        self.time = time
        pos_vec, vel_vec = (
            self.input_coords.si_values()[:3],
            self.input_coords.si_values()[3:],
        )
        time_vel = schwarzschild_utils.time_velocity(pos_vec, vel_vec, M)
        self.initial_vec = np.hstack(
            (self.time.value, pos_vec, time_vel.value, vel_vec)
        )
        self.scr = schwarzschild_radius(M)

    @classmethod
    @u.quantity_input(time=u.s, M=u.kg, a=u.m)
    def from_coords(cls, coords, M, q=None, Q=None, time=0 * u.s, a=0 * u.m):
        """
        Constructor

        Parameters
        ----------
        coords : ~einsteinpy.coordinates.velocity.CartesianDifferential
            Object having both initial positions and velocities of particle in Cartesian Coordinates
        M : ~astropy.units.quantity.Quantity
            Mass of the body
        time : ~astropy.units.quantity.Quantity
            Time of start, defaults to 0 seconds.

        """
        if coords.system == "Spherical":
            return cls(coords, M, time)
        sph_coords = coords.spherical_differential()
        return cls(sph_coords, M, time)

    def f_vec(self, ld, vec):
        vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
        chl = schwarzschild_utils.christoffels(vec[1], vec[2], self.M.value)
        vals[:4] = vec[4:8]
        vals[4] = -2 * chl[0, 0, 1] * vec[4] * vec[5]
        vals[5] = -1 * (
            chl[1, 0, 0] * (vec[4] ** 2)
            + chl[1, 1, 1] * (vec[5] ** 2)
            + chl[1, 2, 2] * (vec[6] ** 2)
            + chl[1, 3, 3] * (vec[7] ** 2)
        )
        vals[6] = -2 * chl[2, 2, 1] * vec[6] * vec[5] - 1 * chl[2, 3, 3] * (vec[7] ** 2)
        vals[7] = -2 * (chl[3, 3, 1] * vec[7] * vec[5] + chl[3, 3, 2] * vec[7] * vec[6])
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
            Starting lambda(proper time), defaults to 0, (lambda ~= t)
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
        ODE = RK45(
            fun=self.f_vec,
            t0=start_lambda,
            y0=self.initial_vec,
            t_bound=end_lambda,
            **OdeMethodKwargs
        )

        vecs = list()
        lambdas = list()
        crossed_event_horizon = False
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
                    SphericalDifferential(
                        v[1] * u.m,
                        v[2] * u.rad,
                        v[3] * u.rad,
                        v[5] * u.m / u.s,
                        v[6] * u.rad / u.s,
                        v[7] * u.rad / u.s,
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
                    SphericalDifferential(
                        v[1] * u.m,
                        v[2] * u.rad,
                        v[3] * u.rad,
                        v[5] * u.m / u.s,
                        v[6] * u.rad / u.s,
                        v[7] * u.rad / u.s,
                    )
                    .cartesian_differential()
                    .si_values()
                )
                yield ODE.t, np.hstack((v[0], si_vals[:3], v[4], si_vals[3:]))
            ODE.step()
            if (not crossed_event_horizon) and (ODE.y[1] <= _scr):
                warnings.warn("particle reached Schwarzchild Radius. ", RuntimeWarning)
                if stop_on_singularity:
                    break
                else:
                    crossed_event_horizon = True
