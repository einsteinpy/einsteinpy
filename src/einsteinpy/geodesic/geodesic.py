import warnings

import numpy as np
from scipy.optimize import fsolve

from .utils import _energy, _julia_solver, _python_solver, _sphToCart


class Geodesic:
    """
    Base Class for defining Geodesics
    Working in Geometrized Units (M-Units), with ``G = c = M = 1.``

    """

    def __init__(
        self,
        position,
        momentum,
        a=0.0,
        end_lambda=50.0,
        step_size=0.0005,
        time_like=True,
        return_cartesian=True,
        julia=True,
    ):
        """
        Constructor

        Parameters
        ----------
        position : array_like
            Length-3 Array, containing the initial 3-Position
        momentum : array_like
            Length-3 Array, containing the initial 3-Momentum
        a : float, optional
            Dimensionless Spin Parameter of the Black Hole
            ``0 <= a <= 1``
            Defaults to ``0.`` (Schwarzschild Black Hole)
        end_lambda : float, optional
            Affine Parameter value, where integration will end
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``50.``
        step_size : float, optional
            Size of each geodesic integration step
            A fixed-step, symplectic VerletLeapfrog integrator is used
            Defaults to ``0.0005``
        time_like : bool, optional
            Determines type of Geodesic
            ``True`` for Time-like geodesics
            ``False`` for Null-like geodesics
            Defaults to ``True``
        return_cartesian : bool, optional
            Whether to return calculated positions in Cartesian Coordinates
            This only affects the coordinates. The momenta dimensionless quantities,
            and are returned in Spherical Polar Coordinates.
            Defaults to ``True``
        julia : bool, optional
            Whether to use the julia backend
            Defaults to ``True``

        """
        self.position = position
        self.momentum = momentum
        self.a = a
        self.end_lambda = end_lambda
        self.step_size = step_size
        self.kind = "Time-like" if time_like else "Null-like"
        self.coords = "Cartesian" if return_cartesian else " Spherical Polar"
        self.backend = "Julia" if julia else "Python"

        self._trajectory = self.calculate_trajectory()

    def __repr__(self):
        return f"Geodesic Object:\n\
            Type = ({self.kind}),\n\
            Position = ({self.position}),\n\
            Momentum = ({self.momentum}),\n\
            Spin Parameter = ({self.a})\n\
            Solver details = (\n\
                Backend = ({self.backend})\n\
                Step-size = ({self.step_size}),\n\
                End-Lambda = ({self.end_lambda})\n\
                Trajectory = (\n\
                    {self.trajectory}\n\
                ),\n\
                Output Position Coordinate System = ({self.coords})\n\
            )"

    def __str__(self):
        return self.__repr__()

    @property
    def trajectory(self):
        """
        Returns the trajectory of the test particle

        """
        return self._trajectory

    def calculate_trajectory(self):
        """
        Calculate trajectory in spacetime, according to Geodesic Equations

        Returns
        -------
        ~numpy.ndarray
            N-element numpy array, containing affine parameter 
            values, where the integration was performed
        ~numpy.ndarray
            Shape-(N, 6) numpy array, containing [x1, x2, x3, p_r, p_theta, p_phi] for each Lambda

        """
        mu = 1.0 if self.kind == "Time-like" else 0.0
        q, p = self.position, self.momentum
        a = self.a

        # Getting Energy value, after solving guu.pd.pd = -mu ** 2, where,
        # 'u' denotes contravariant index and 'd' denotes covariant index
        E = fsolve(_energy, 1.0, args=(q, p, a, mu))[-1]

        params = [a, E, mu]

        if self.backend == "Python":
            warnings.warn(
                """
                Using Python backend to solve the system. This backend is currently in beta and the solver 
                may not be stable for certain sets of conditions, e.g. long simulations (`end_lambda > 50.`) 
                or high initial radial distances (`position[0] > ~5.`). In these cases or if the output does not seem accurate,
                it is highly recommended to switch to the Julia backend, by setting `julia=True`, in the constructor call.
                """,
                RuntimeWarning,
            )
            lambdas, vecs = _python_solver(q, p, params, end_lambda, step_size)

        else:
            lambdas, vecs = _julia_solver(q, p, params, end_lambda, step_size)

        if self.coords == "Cartesian":
            xc = list()
            yc = list()
            zc = list()

            # Converting to Cartesian from Spherical Polar Coordinates
            # Note that momenta cannot be converted correctly, this way,
            # due to ambiguities in the signs of v_r and v_th (velocities)
            cart_vecs = list()
            for y in vecs:
                r, th, phi = y[0], y[1], y[2]
                cart_coords = _sphToCart(r, th, phi)
                cart_vecs.append(np.hstack((cart_coords, y[3:])))

            return lambdas, np.array(cart_vecs)

        return lambdas, vecs
