import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.metric import BaseMetric
from einsteinpy.units import primitive
from einsteinpy.utils import CoordinateError

_c = constant.c.value
_G = constant.G.value
_eps0 = constant.eps0.value


class ReissnerNordstrom(BaseMetric):
    """
    Class for defining Reissner-Nordstrom Geometry
    using the (+,-,-,-) signature

    """

    @u.quantity_input(M=u.kg, Q=u.C)
    def __init__(self, coords, M, Q):
        """
        Constructor

        Parameters
        ----------
        coords : ~einsteinpy.coordinates.differential.*
            Coordinate system, in which Metric is to be represented
        M : ~astropy.units.quantity.Quantity
            Mass of gravitating body, e.g. Black Hole
        Q : ~astropy.units.quantity.Quantity
            Charge on gravitating body, e.g. Black Hole

        """
        super().__init__(
            coords=coords,
            M=M,
            Q=Q,
            name="ReissnerNordstrom Metric",
            metric_cov=self.metric_covariant,
            christoffels=self._christoffels,
            f_vec=self._f_vec,
        )

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Reissner-Nordstrom Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Reissner-Nordstrom Metric Tensor in chosen Coordinates
            Numpy array of shape (4,4) with signature (+,-,-,-)

        """
        if self.coords.system == "Spherical":
            return self._g_cov_s(x_vec)

        raise CoordinateError(
            "Reissner-Nordstrom Metric is available only in spherical polar coordinates."
        )

    def _g_cov_s(self, x_vec):
        """
        Returns Covariant Reissner-Nordstrom Metric Tensor \
        in spherical coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Reissner-Nordstrom Metric Tensor
            Numpy array of shape (4,4) with convention (+,-,-,-)

        """
        r, th = x_vec[1], x_vec[2]

        # define the constants used in the metric
        r_s = self.sch_rad
        rho = _G / (_eps0 * _c**4)
        (_Q,) = primitive(self.Q)

        # as the metric is diagonal, initialise with a zero metric
        g_cov = np.zeros(shape=(4, 4), dtype=float)

        # then define the diagonal terms
        Arn = 1 - r_s / r + (rho * _Q**2) / (r**2)
        g_cov[0, 0] = Arn * _c**2
        g_cov[1, 1] = -1.0 / Arn
        g_cov[2, 2] = -(r**2)
        g_cov[3, 3] = -((r * np.sin(th)) ** 2)

        return g_cov

    def _christoffels(self, x_vec):
        """
        Returns Christoffel Symbols for Reissner-Nordstrom Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Reissner-Nordstrom Metric \
            in chosen Coordinates
            Numpy array of shape (4,4,4)

        Raises
        ------
        CoordinateError
            Raised, if the Christoffel Symbols are not \
            available in the supplied Coordinate System

        """
        if self.coords.system == "Spherical":
            return self._ch_sym_s(x_vec)

        raise CoordinateError(
            "Christoffel Symbols for Reissner-Nordstrom Metric are available only in Spherical Polar Coordinates."
        )

    def _ch_sym_s(self, x_vec):
        """
        Returns the Christoffel Symbols for Reissner-Nordstrom metric
        in spherical coordinates, with (u,l,l) indices

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Reissner-Nordstrom Metric \
            in spherical coordinates, with (u,l,l) indices
            Numpy array of shape (4,4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s = self.sch_rad
        rho = _G / (_eps0 * _c**4)
        (_Q,) = primitive(self.Q)
        Arn = 1 - r_s / r + (rho * _Q**2) / (r**2)

        # set all elements to zero then reset the non-zero terms
        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        # \Gamma(t,t,r)
        chl[0, 0, 1] = chl[0, 1, 0] = (r_s * r - 2 * rho * _Q**2) / (2 * r**3 * Arn)

        # \Gamma(r,t,t)
        chl[1, 0, 0] = Arn * (_c**2) * (r_s * r - 2 * rho * _Q**2) / (2 * r**3)

        # \Gamma(r,r,r)
        chl[1, 1, 1] = -(r_s * r - 2 * rho * _Q**2) / (2 * r**3 * Arn)

        # \Gamma(r,th,th)
        chl[1, 2, 2] = -r * Arn

        # \Gamma(r,phi,phi)
        chl[1, 3, 3] = -r * Arn * (np.sin(th) ** 2)

        # \Gamma(th,th,r)
        chl[2, 1, 2] = chl[2, 2, 1] = 1 / r

        # \Gamma(th,phi,phi)
        chl[2, 3, 3] = -np.sin(th) * np.cos(th)

        # \Gamma(phi,phi,r)
        chl[3, 3, 1] = chl[3, 1, 3] = 1 / r

        # \Gamma(phi,phi,th)
        chl[3, 3, 2] = chl[3, 2, 3] = 1 / np.tan(th)

        return chl

    def _f_vec(self, lambda_, vec):
        """
        Returns f_vec for Reissner-Nordstrom Metric in chosen coordinates
        To be used for solving Geodesics ODE

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
            f_vec for Kerr-Newman Metric in chosen coordinates
            Numpy array of shape (8)

        Raises
        ------
        CoordinateError
            Raised, if ``f_vec`` is not available in \
            the supplied Coordinate System

        """
        if self.coords.system == "Spherical":
            return self._f_vec_s(lambda_, vec)

        raise CoordinateError(
            "'f_vec' for Schwarzschild Metric is available only in Spherical Polar Coordinates."
        )

    def _f_vec_s(self, lambda_, vec):
        """
        Returns f_vec for the Reissner-Nordstrom Metric \
        in spherical Coordinates
        To be used for solving Geodesics ODE

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
            f_vec for Reissner-NOrdstrom Metric in spherical Coordinates
            Numpy array of shape (8)

        """
        chl = self.christoffels(vec[:4])

        vals = np.zeros(shape=vec.shape, dtype=vec.dtype)

        vals[:4] = vec[4:]

        vals[:4] = vec[4:]
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
