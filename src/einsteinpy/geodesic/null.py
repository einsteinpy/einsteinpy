import warnings

import numpy as np
from scipy.integrate import odeint

from einsteinpy.coordinates.utils import spherical_to_cartesian_fast
from einsteinpy.geodesic.utils import _f_vec, _g_dd, _v0


class Nulllike:
    """
    Class for calculating Null-like Geodesics in Kerr Spacetime, \
    in Boyer-Lindquist Coordinates and M-Units (G = c = M = 1)

    """

    def __init__(
        self,
        position,
        velocity,
        a,
        end_lambda=200.0,
        max_steps=200,
        return_cartesian=True,
    ):
        """
        Parameters
        ----------
        position : array_like
            4-Position of the test particle, in Boyer-Lindquist \
            Coordinates, [t, r, theta, phi]
        velocity : array_like
            Proper velocities, of the test particle i.e. velocities, \
            with respect to affine parameter, in Boyer-Lindquist \
            Coordinates, [v_r, v_theta, v_phi]
        a : float
            Spin Parameter, 0 <= a <= 1
        end_lambda : float, optional
            Value of affine parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``200.``
        max_steps : int, optional
            Maximum number of steps, that the integrator takes
            Defaults to ``200``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Raises
        ------
        ValueError
            If a is not between 0 and 1

        """
        if a < 0 or a > 1:
            raise ValueError("Spin Parameter, a, should be between 0 and 1.")

        self.position = position
        self.velocity = velocity
        self.a = a

        self._state = self._calculate_state()
        self._trajectory = self.calculate_trajectory(
            end_lambda=end_lambda,
            max_steps=max_steps,
            return_cartesian=return_cartesian,
        )

    def __repr__(self):
        return f"NullGeodesic Object:\n\
            Black Hole Spin = ({self.a})\n\
            Initial State = ({self.state}),\n\
            Trajectory = ({self.trajectory})"

    def __str__(self):
        return f"NullGeodesic Object:\n\
            Black Hole Spin = ({self.a})\n\
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
        Returns the trajectory of the Geodesic

        """
        return self._trajectory

    def _calculate_state(self):
        """
        Prepares and returns the Initial State Vector of the massless test particle

        Returns
        -------
        state : ~numpy.ndarray
            Initial State Vector of the massless test particle
            Length-8 Array

        """
        pos, vel = self.position, self.velocity

        g = _g_dd(pos[1], pos[2], self.a)
        v1, v2, v3 = vel

        v_t = _v0(g, v1, v2, v3)
        state = np.hstack((pos, v_t, vel))

        return state

    def calculate_trajectory(
        self, end_lambda=200.0, max_steps=200, return_cartesian=True
    ):
        """
        Solves the Geodesic Equation using `scipy.integrate.odeint`

        Parameters
        ----------
        end_lambda : float, optional
            Value of affine parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``200.``
        max_steps : int, optional
            Maximum number of steps, that the integrator takes
            Defaults to ``200``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray
            N-element numpy array containing Lambda, where the geodesic equations were evaluated
            Here, `N == max_steps`.
        ~numpy.ndarray
            (N,8) shape numpy array containing [x0, x1, x2, x3, v0, v1, v2, v3] for each Lambda
            Here, `N == max_steps`.

        """
        y = self.state
        lambdas = np.linspace(0, end_lambda, max_steps)

        sol = odeint(_f_vec, y, lambdas, args=(self.a,), tfirst=True)

        if return_cartesian:
            cart_vecs = list()
            for vs in sol:
                vc = spherical_to_cartesian_fast(
                    vs[0],
                    vs[1],
                    vs[2],
                    vs[3],
                    vs[5],
                    vs[6],
                    vs[7],
                    velocities_provided=True,
                )
                cart_vecs.append(np.hstack((vc[:4], vs[4], vc[4:])))

            return lambdas, np.array(cart_vecs)

        return lambdas, sol
