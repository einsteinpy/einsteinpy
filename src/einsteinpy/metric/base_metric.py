"""Docstring for base_metric.py module

This module defines the ``BaseMetric`` class, which is the base class for all \
metrics in EinsteinPy. This class contains several utilities, that are \
used in ``einsteinpy.metric`` to define classes for vacuum solutions. Users are \
advised to inherit this class, while defining their own metric classes. \
Two parameters to note are briefly described below:

- ``metric_cov``: User-supplied function, defining the covariant form of the metric tensor. Users need \
    to supply just this to completely determine the metric tensor, as the contravariant form \
    is calculated and accessed through a predefined method, ``metric_contravariant()``.
- ``perturbation``: User-supplied function, defining a perturbation to the metric. Currently, no checks are \
    performed to ascertain the physicality of the resulting perturbed metric. Read the documentation \
    on ``metric_covariant()`` below, to learn more.

Also, note that, users should call ``metric_covariant`` to access the perturbed, covariant form of the metric. \
For unperturbed underlying metric, users should call ``metric_cov``, which returns the metric, that they had \
supplied.

"""
import warnings

import numpy as np

from einsteinpy import constant

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


class BaseMetric:
    """
    For defining a general Metric

    """

    def __init__(
        self,
        coords,
        M,
        a=0,
        Q=0,
        name="Base Metric",
        metric_cov=None,
        christoffels=None,
        f_vec=None,
        perturbation=None,
    ):
        """
        Constructor

        Parameters
        ----------
        coords : string
            Coordinate system, in which Metric is to be represented
            "S" - Schwarzschild: Only applicable to Schwarzschild solutions
            "BL" - Boyer-Lindquist: Applicable to Kerr-Newman solutions
            "KS" - Kerr-Schild: Useful for adding perturbations to Kerr-Newman solutions
        M : float
            Mass of gravitating body, e.g. Black Hole
        a : float
            Spin Parameter
            Defaults to ``0``
        Q : float
            Charge on gravitating body, e.g. Black Hole
            Defaults to ``0``
        name : str, optional
            Name of the Metric Tensor. Defaults to ``"Base Metric"``
        metric_cov : callable, optional
            Function, defining Covariant Metric Tensor
            It should return a real-valued tensor (2D Array), at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        christoffels : callable, optional
            Function, defining Christoffel Symbols
            It should return a real-valued (4,4,4) array, at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        f_vec : callable, optional
            Function, defining RHS of Geodesic Equation
            It should return a real-valued (8) vector, at supplied coordinates
            Defaults to ``None``
            Consult pre-defined classes for function definition
        perturbation : callable, optional
            Function, defining Covariant Perturbation Tensor
            It should return a real-valued tensor (2D Array), at supplied coordinates
            Defaults to ``None``
            Function definition similar to ``metric_cov``

        """
        self.name = name
        self.coords = coords
        self.M = M
        self.a = a
        self.Q = Q
        self.metric_cov = metric_cov
        self.christoffels = christoffels
        self.f_vec = f_vec
        self.perturbation = perturbation
        # Need physical checks for `perturbation`
        # Expert opinion needed - Gauge Fixing

        self.sch_rad = self.schwarzschild_radius(M)

    def __str__(self):
        return f"(\nName: ({self.name}),\
            \nCoordinates: ({self.coords}),\
            \nMass: ({self.M}),\
            \nSpin parameter: ({self.a}),\
            \nCharge: ({self.Q}),\
            \nSchwarzschild Radius: ({self.sch_rad})\n)"

    def __repr__(self):
        return f"(\nName: ({self.name}),\
            \nCoordinates: ({self.coords}),\
            \nMass: ({self.M}),\
            \nSpin parameter: ({self.a}),\
            \nCharge: ({self.Q}),\
            \nSchwarzschild Radius: ({self.sch_rad})\n)"

    @staticmethod
    def sigma(r, theta, M, a):
        """
        Returns the value of ``r**2 + alpha**2 * cos(theta)**2``
        Specific to Boyer-Lindquist coordinates
        Applies to Kerr Geometry

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter

        Returns
        -------
        float
            The value of ``r**2 + alpha**2 * cos(theta)**2``

        """
        alpha = BaseMetric.alpha(M, a)
        sigma = (r ** 2) + ((alpha * np.cos(theta)) ** 2)
        return sigma

    @staticmethod
    def delta(r, M, a, Q=0):
        """
        Returns the value of ``r**2 - r_s * r + alpha**2 + r_Q**2``
        Specific to Boyer-Lindquist coordinates
        Applies to Kerr & Kerr-Newman Geometries

        Parameters
        ----------
        r : float
            r-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float, optional
            Charge on gravitating body
            Defaults to ``0`` (for Kerr Geometry)

        Returns
        -------
        float
            The value of ``r**2 - r_s * r + alpha**2 + r_Q**2``

        """
        r_s = 2 * M * _G / _c ** 2
        alpha = BaseMetric.alpha(M, a)
        # Square of Geometrized Charge
        r_Q2 = (Q ** 2) * _G * _Cc / _c ** 4
        delta = (r ** 2) - (r_s * r) + (alpha ** 2) + r_Q2
        return delta

    @staticmethod
    def rho(r, theta, M, a):
        """
        Returns the value of ``sqrt(r**2 + alpha**2 * cos(theta)**2) == sqrt(sigma)``
        Specific to Boyer-Lindquist coordinates
        Applies to Kerr-Newman Geometry

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter

        Returns
        -------
        float
            The value of ``sqrt(r**2 + alpha**2 * cos(theta)**2) == sqrt(sigma)``

        """
        return np.sqrt(BaseMetric.sigma(r, theta, M, a))

    @staticmethod
    def r_ks(x, y, z, a):
        """
        Returns the value of r, after solving \
        ``(x**2 + y**2) / (r**2 + alpha**2) + z**2 / r**2 = 1``
        'r' is not the Radius Coordinate of Spherical Polar or BL Coordinates
        Specific to Cartesian form of Kerr-Schild Coordinates
        """
        raise NotImplementedError

    @staticmethod
    def schwarzschild_radius(M):
        """
        Returns Schwarzschild Radius

        Parameters
        ----------
        M : float
            Mass of gravitating body

        Returns
        -------
        r : float
            Schwarzschild Radius for a given mass
        """
        return 2 * _G * M / _c ** 2

    @staticmethod
    def alpha(M, a):
        """
        Returns Rotational Length Parameter (alpha) that is used \
        in the Metric. Following equations are relevant:
        alpha = J / Mc
        a = Jc / GM^2
        alpha = GMa / c^2
        where, 'a' is the dimensionless Spin Parameter (0 <= a <= 1)

        Parameters
        ----------
        M : float
            Mass of gravitating body
        a : float
            Number between 0 and 1

        Returns
        -------
        float
            Rotational Length Parameter

        Raises
        ------
        ValueError
            If a is not between 0 and 1

        """
        if a < 0 or a > 1:
            raise ValueError("Spin Parameter, a, should be between 0 and 1.")
        half_rs = _G * M / _c ** 2

        return a * half_rs

    def singularities(self):
        """
        Returns the Singularities of the Metric
        Depends on the choice of Coordinate Systems
        Applies to Kerr and Kerr-Newman Geometries

        Returns
        -------
        dict
            Dictionary of singularities in the geometry
            ``{
            "inner_ergosphere": function(theta),
            "inner_horizon": float,
            "outer_horizon": float,
            "outer_ergosphere": function(theta)
            }``

        """
        coords, M, a, Q = self.coords, self.M, self.a, self.Q

        r_s = 2 * M * _G / _c ** 2
        alpha = BaseMetric.alpha(M, a)
        # Square of Geometrized Charge
        r_Q2 = (Q ** 2) * _G * _Cc / _c ** 4

        def _in_ergo(theta):
            return (
                r_s
                - np.sqrt((r_s ** 2) - (4 * (alpha * np.cos(theta)) ** 2) - (4 * r_Q2))
            ) / 2

        def _out_ergo(theta):
            return (
                r_s
                + np.sqrt((r_s ** 2) - (4 * (alpha * np.cos(theta)) ** 2) - (4 * r_Q2))
            ) / 2

        if coords == "S":  # Schwarzschild Geometry
            return {
                "inner_ergosphere": 0,
                "inner_horizon": 0,
                "outer_horizon": r_s,
                "outer_ergosphere": r_s,
            }

        elif coords == "BL":  # Kerr & Kerr-Newman Geometries
            return {
                "inner_ergosphere": _in_ergo,
                "inner_horizon": (
                    r_s - np.sqrt((r_s ** 2) - (4 * alpha ** 2) - (4 * r_Q2))
                )
                / 2,
                "outer_horizon": (
                    r_s + np.sqrt((r_s ** 2) - (4 * alpha ** 2) - (4 * r_Q2))
                )
                / 2,
                "outer_ergosphere": _out_ergo,
            }

        elif coords == "KS":  # Kerr & Kerr-Newman Geometries
            # To be filled in, after refactoring `coordinates`
            raise NotImplementedError

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Metric Tensor
        Adds Kerr-Schild (Linear) Perturbation to metric, \
        if ``perturbation`` is defined in Metric object
        Currently, this does not consider Gauge Fixing or \
        any physical checks for the ``perturbation`` matrix. \
        Please exercise caution while using ``perturbation``.

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Metric Tensor
            Numpy array of shape (4,4)

        """
        g_cov = self.metric_cov(x_vec)

        if self.perturbation:
            p_cov = self.perturbation(x_vec)
            return g_cov + p_cov

        return g_cov

    def metric_contravariant(self, x_vec):
        """
        Returns Contravariant Metric Tensor
        Adds Kerr-Schild (Linear) Perturbation to metric, \
        if ``perturbation`` is not None in Metric object
        Currently, this does not consider Gauge Fixing or \
        any physical checks for the ``perturbation`` matrix. \
        Please exercise caution while using ``perturbation``.

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Contravariant Metric Tensor

        """
        g_cov = self.metric_covariant(x_vec)
        g_contra = np.linalg.inv(g_cov)

        return g_contra

    def em_potential_covariant(self, r, theta, M, a, Q):
        """
        Returns Covariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body

        Returns
        -------
        ~numpy.ndarray
            Covariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        alpha = BaseMetric.alpha(M, a)
        # Geometrized Charge
        r_Q = np.sqrt((Q ** 2 * _G * _Cc) / _c ** 4)
        rho2 = BaseMetric.rho(r, theta, M, a) ** 2

        A = np.zeros((4,), dtype=float)
        A[0] = r * r_Q / rho2
        A[3] = -r * alpha * r_Q * np.sin(theta) ** 2 / rho2

        return A

    def em_potential_contravariant(self, r, theta, M, a, Q):
        """
        Returns Contravariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body

        Returns
        -------
        ~numpy.ndarray
            Contravariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        A_cov = self.em_potential_covariant(r, theta, M, a, Q)
        x_vec = np.array(
            [0.0, r, theta, 0.0], dtype=float
        )  # t & phi have no bearing on Metric
        g_contra = self.metric_contravariant(x_vec=x_vec)

        return g_contra @ A_cov

    def em_tensor_covariant(self, r, theta, M, a, Q):
        """
        Returns Covariant Electromagnetic Tensor
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body

        Returns
        -------
        ~numpy.ndarray
            Covariant Electromagnetic Tensor
            Numpy array of shape (4, 4)

        """
        alpha = BaseMetric.alpha(M, a)
        r_Q = np.sqrt((Q ** 2 * _G * _Cc) / _c ** 4)
        rho2 = BaseMetric.rho(r, theta, M, a) ** 2
        # Partial derivatives of rho2
        drho2_dr = 2 * r
        drho2_dtheta = -(alpha ** 2 * np.sin(2 * theta))

        F = np.zeros((4, 4), dtype=float)

        F[0, 1] = -(r_Q * (rho2 - drho2_dr * r)) / (rho2 ** 2)
        F[1, 0] = -F[0, 1]
        F[0, 2] = (r * r_Q * drho2_dtheta) / (rho2 ** 2)
        F[2, 0] = -F[0, 2]
        F[1, 3] = (
            (1 / rho2 ** 2) * (alpha * r_Q * np.sin(theta) ** 2) * (rho2 - 2 * r ** 2)
        )
        F[3, 1] = -F[1, 3]
        F[2, 3] = (
            (1 / rho2 ** 2)
            * (alpha * r_Q * r * np.sin(2 * theta))
            * (rho2 + (alpha * np.sin(theta)) ** 2)
        )
        F[3, 2] = -F[2, 3]

        return F

    def em_tensor_contravariant(self, r, theta, M, a, Q):
        """
        Returns Contravariant Electromagnetic Tensor
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        r : float
            r-component of 4-Position
        theta : float
            theta-component of 4-Position
        M : float
            Mass of gravitating body
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body

        Returns
        -------
        ~numpy.ndarray
            Contravariant Electromagnetic Tensor
            Numpy array of shape (4, 4)

        """
        F_cov = self.em_tensor_covariant(r, theta, M, a, Q)
        x_vec = np.array([0, r, theta, 0], dtype=float)
        g_contra = self.metric_contravariant(x_vec=x_vec)

        F_contra = g_contra @ F_cov @ g_contra

        return F_contra

    # Deprecated function
    def calculate_trajectory(
        self,
        start_lambda=0.0,
        end_lambda=10.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        return_cartesian=False,
    ):
        """
        Deprecated in 0.4.0.
        Please use ``einsteinpy.Geodesic``.

        Calculate trajectory in manifold according to geodesic equation

        """
        warnings.warn(
            "calculate_trajectory() \
            has been deprecated in Version 0.4.0 \
            Please use einsteinpy.Geodesic.calculate_trajectory()!",
            DeprecationWarning,
        )
