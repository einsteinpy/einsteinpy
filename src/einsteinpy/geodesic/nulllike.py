import warnings

import numpy as np

from einsteinpy.coordinates.utils import bl_to_cartesian_fast
from einsteinpy.geodesic.utils import _calc_state, _christoffels_M, _f_geod_M
from einsteinpy.integrators import RK45


class NullGeodesic:
    """
    Class for defining Null-like Geodesics

    """

    # - ????? Too many parameters -
    # Can remove E
    # Consolidate alpha, beta, rcam & inclination
    # Consolidate end_lambda, step_size, return_cartesian as ODEKwargs
    def __init__(
        self,
        alpha,
        beta,
        rcam=1e4,
        inclination=np.pi / 2,
        a=0.0,
        E=1.0,
        end_lambda=10.0,
        step_size=1e-3,
        return_cartesian=True,
    ):
        """
        Parameters
        ----------
        # - Docstring ?????
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
        step_size : float, optional
            Size of each geodesic integration step
            Defaults to ``1e-3``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        """
        self.alpha = alpha
        self.beta = beta
        self.rcam = rcam
        self.inclination = inclination
        self.a = a
        self.E = E  # - ?????

        # self._sch_rad = 2
        self._state = self._calculate_state()
        self._trajectory = self.calculate_trajectory(
            end_lambda=end_lambda,
            OdeMethodKwargs={"stepsize": step_size},
            return_cartesian=return_cartesian,
        )[1]

    def __repr__(self):
        return f"NullGeodesic Object:\n\
            Initial Conditions = ({self.alpha, self.beta, self.rcam, self.inclination}),\n\
            Initial State = ({self.state}),\n\
            Trajectory = ({self.trajectory})"

    def __str__(self):
        return f"NullGeodesic Object:\n\
            Initial Conditions = ({self.alpha, self.beta, self.rcam, self.inclination}),\n\
            Initial State = ({self.state}),\n\
            Trajectory = ({self.trajectory})"

    @property
    def state(self):
        """
        Returns the Initial State Vector of the Geodesic

        """
        return self._state

    @property
    def trajectory(self):
        """
        Returns the "Trajectory" of the Geodesic

        """
        return self._trajectory

    def _calculate_state(self):
        """
        Prepares and returns the Initial State Vector of the massless test particle

        Source: RAPTOR (?????)

        Returns
        -------
        state : ~numpy.ndarray
            Initial State Vector of the massless test particle
            Length-8 Array

        """
        a = self.a
        E = self.E
        r = self.rcam
        alpha, beta = self.alpha, self.beta
        theta = self.inclination

        return _calc_state(alpha, beta, r, theta, a, E)

    # Move to Kerr - ?????
    # Not until units + coordinate switching is figured out
    def _ch_sym_M(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr Metric \
        in Boyer-Lindquist Coordinates, in M-Units

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr Metric \
            in Boyer-Lindquist Coordinates \
            in M-Units
            Numpy array of shape (4,4,4)

        """
        r, th = x_vec[1], x_vec[2]

        return _christoffels_M(self.a, r, th)

    # Move to Kerr - ?????
    # Not until units + coordinate switching is figured out
    def _f_vec_M(self, lambda_, vec):
        """
        Returns f_vec for Kerr Metric \
        in Boyer-Lindquist Coordinates \
        in M-Units (G = c = M = 1)

        To be used for solving Geodesics ODE

        Source: RAPTOR (?????)

        Parameters
        ----------
        lambda_ : float
            Parameterizes current integration step
            Used by ODE Solver

        vec : array_like
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.ndarray
            f_vec for Kerr Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)

        """
        chl = self._ch_sym_M(vec[:4])
        vals = np.zeros(shape=vec.shape, dtype=vec.dtype)

        return _f_geod_M(chl, vals, vec)

    def calculate_trajectory(
        self,
        end_lambda=10.0,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=True,
    ):
        """
        Calculate trajectory in spacetime, according to Geodesic Equations

        Parameters
        ----------
        end_lambda : float, optional
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``10.0``
        OdeMethodKwargs : dict, optional
            Kwargs to be supplied to the ODESolver
            Dictionary with key 'stepsize' along with a float value is expected
            Defaults to ``{'stepsize': 1e-3}``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray
            N-element numpy array containing Lambda, where the geodesic equations were evaluated
        ~numpy.ndarray
            (n,8) shape numpy array containing [x0, x1, x2, x3, v0, v1, v2, v3] for each Lambda

        """
        ODE = RK45(
            fun=self._f_vec_M,
            t0=0.0,
            y0=self.state,
            t_bound=end_lambda,
            **OdeMethodKwargs,
        )

        a = self.a
        r = self.rcam

        # Termination conditions
        cutoff_outer = r * 1.01
        cutoff_inner = (1.0 + np.sqrt(1.0 - a ** 2)) * 1.1

        vecs = list()
        lambdas = list()

        while ODE.t < end_lambda:
            vecs.append(ODE.y)
            lambdas.append(ODE.t)

            r_curr = ODE.y[1]

            # Checking termination conditions
            if r_curr < cutoff_inner or r_curr > cutoff_outer:
                warnings.warn("Light Ray has reached cut-off bounds.", RuntimeWarning)
                break

            ODE.step()

        vecs, lambdas = np.array(vecs), np.array(lambdas)

        if return_cartesian:
            cart_vecs = list()
            for v in vecs:
                vals = bl_to_cartesian_fast(
                    v[0],
                    v[1],
                    v[2],
                    v[3],
                    a,
                    v[5],
                    v[6],
                    v[7],
                    velocities_provided=True,
                )
                cart_vecs.append(np.hstack((vals[:4], v[4], vals[4:])))

            return lambdas, np.array(cart_vecs)

        return lambdas, vecs


class NullBundle:
    """
    Class for generating a photon sheet and performing
    geodesic integration for Radiative Transfer applications

    """
