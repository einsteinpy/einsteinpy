import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.metric import BaseMetric
from einsteinpy.utils import CoordinateError

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


class KerrNewman(BaseMetric):
    """
    Class for defining Kerr-Newman Goemetry

    """

    @u.quantity_input(M=u.kg, a=u.one, Q=u.C, q=u.C / u.kg)
    def __init__(self, coords, M, a, Q, q=0.0 * u.C / u.kg):
        """
        Constructor

        Parameters
        ----------
        coords : ~einsteinpy.coordinates.differential.*
            Coordinate system, in which Metric is to be represented
        M : ~astropy.units.quantity.Quantity
            Mass of gravitating body, e.g. Black Hole
        a : ~astropy.units.quantity.Quantity
            Spin Parameter
        Q : ~astropy.units.quantity.Quantity
            Charge on gravitating body, e.g. Black Hole
        q : ~astropy.units.quantity.Quantity, optional
            Charge, per unit mass, of the test particle
            Defaults to ``0 C / kg``

        """
        super().__init__(
            coords=coords,
            M=M,
            a=a,
            Q=Q,
            name="Kerr-Newman Metric",
            metric_cov=self.metric_covariant,
            christoffels=self._christoffels,
            f_vec=self._f_vec,
        )
        self.q = q
        # Precomputed list of tuples, containing indices \
        # of non-zero Christoffel Symbols for Kerr-Newman Metric \
        # in Boyer-Lindquist Coordinates
        self._nonzero_christoffels_list_bl = [
            (0, 0, 1),
            (0, 0, 2),
            (0, 1, 3),
            (0, 2, 3),
            (0, 1, 0),
            (0, 2, 0),
            (0, 3, 1),
            (0, 3, 2),
            (1, 0, 0),
            (1, 1, 1),
            (1, 2, 2),
            (1, 3, 3),
            (2, 0, 0),
            (2, 1, 1),
            (2, 2, 2),
            (2, 3, 3),
            (1, 0, 3),
            (1, 1, 2),
            (2, 0, 3),
            (2, 1, 2),
            (1, 2, 1),
            (1, 3, 0),
            (2, 2, 1),
            (2, 3, 0),
            (3, 0, 1),
            (3, 0, 2),
            (3, 1, 0),
            (3, 1, 3),
            (3, 2, 0),
            (3, 2, 3),
            (3, 3, 1),
            (3, 3, 2),
        ]

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr-Newman Metric Tensor in chosen Coordinates
            Numpy array of shape (4,4)

        Raises
        ------
        CoordinateError
            Raised, if the metric is not available in \
            the supplied Coordinate System

        """
        if self.coords.system == "BoyerLindquist":
            return self._g_cov_bl(x_vec)

        raise CoordinateError(
            "Kerr-Newman Metric is available only in Boyer-Lindquist Coordinates."
        )

    def _g_cov_bl(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr-Newman Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        M, a = self.M.value, self.a.value
        alpha = super().alpha(M, a)
        rho2, dl = super().rho(r, th, M, a) ** 2, super().delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = (_c**2) * ((dl - ((alpha * np.sin(th)) ** 2)) / (rho2))
        g_cov_bl[1, 1] = -rho2 / dl
        g_cov_bl[2, 2] = -rho2
        g_cov_bl[3, 3] = -(
            (np.sin(th) ** 2)
            * (((r**2 + alpha**2) ** 2 - dl * (alpha * np.sin(th)) ** 2) / rho2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = _c * (
            (-alpha * (np.sin(th) ** 2) * (dl - (r**2) - (alpha**2))) / rho2
        )

        return g_cov_bl

    def _dg_dx_bl(self, x_vec):
        """
        Returns derivative of each Kerr-Newman Metric component \
        w.r.t. coordinates in Boyer-Lindquist Coordinate System

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        dgdx : ~numpy.ndarray
            Array, containing derivative of each Kerr-Newman \
            Metric component w.r.t. coordinates \
            in Boyer-Lindquist Coordinate System
            Numpy array of shape (4,4,4)
            dgdx[0], dgdx[1], dgdx[2] & dgdx[3] contain \
            derivatives of metric w.r.t. t, r, theta & phi respectively

        """
        r, th = x_vec[1], x_vec[2]
        M, a = self.M.value, self.a.value
        alpha = super().alpha(M, a)
        rho2, dl = super().rho(r, th, M, a) ** 2, super().delta(r, M, a)

        dgdx = np.zeros(shape=(4, 4, 4), dtype=float)

        # Metric is invariant on t & phi
        # Differentiation of metric wrt r
        def due_to_r():
            nonlocal dgdx
            drh2dr = 2 * r
            dddr = 2 * r - self.sch_rad
            dgdx[1, 0, 0] = (
                (_c**2)
                * (dddr * rho2 - drh2dr * (dl - (alpha * np.sin(th)) ** 2))
                / (rho2**2)
            )
            dgdx[1, 1, 1] = (-1 / (dl**2)) * (drh2dr * dl - dddr * rho2)
            dgdx[1, 2, 2] = -drh2dr
            dgdx[1, 3, 3] = ((np.sin(th) / rho2) ** 2) * (
                (
                    (
                        ((alpha * np.sin(th)) ** 2) * dddr
                        - 4 * (r**3)
                        - 4 * (r * (alpha**2))
                    )
                    * rho2
                )
                - (
                    drh2dr
                    * (((alpha * np.sin(th)) ** 2) * dl - ((r**2 + alpha**2) ** 2))
                )
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = (
                _c * (-alpha) * (np.sin(th) ** 2) / (rho2**2)
            ) * ((dddr - 2 * r) * rho2 - drh2dr * (dl - r**2 - alpha**2))

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            drh2dth = -(alpha**2) * np.sin(2 * th)
            dgdx[2, 0, 0] = (-((_c / rho2) ** 2)) * (
                (drh2dth * (dl - ((alpha * np.sin(th)) ** 2)))
                + ((alpha**2) * rho2 * np.sin(2 * th))
            )
            dgdx[2, 1, 1] = -drh2dth / dl
            dgdx[2, 2, 2] = -drh2dth
            dgdx[2, 3, 3] = (1 / (rho2**2)) * (
                (dl * (alpha * np.sin(th)) ** 2)
                * (2 * rho2 * np.sin(2 * th) - drh2dth * (np.sin(th)) ** 2)
                - (
                    ((r**2 + alpha**2) ** 2)
                    * (rho2 * np.sin(2 * th) - drh2dth * (np.sin(th)) ** 2)
                )
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = (
                (-alpha * _c * (dl - r**2 - alpha**2)) / (rho2**2)
            ) * ((np.sin(2 * th) * rho2) - (drh2dth * (np.sin(th) ** 2)))

        due_to_r()
        due_to_theta()

        return dgdx

    def _christoffels(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr-Newman \
            Metric in chosen Coordinates
            Numpy array of shape (4,4,4)

        Raises
        ------
        CoordinateError
            Raised, if the Christoffel Symbols are not \
            available in the supplied Coordinate System

        """
        if self.coords.system == "BoyerLindquist":
            return self._ch_sym_bl(x_vec)

        raise CoordinateError(
            "Christoffel Symbols for Kerr-Newman Metric are available only in Boyer-Lindquist Coordinates."
        )

    def _ch_sym_bl(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr-Newman Metric \
            in Boyer-Lindquist Coordinates
            Numpy array of shape (4,4,4)

        """
        g_contra = self.metric_contravariant(x_vec)
        dgdx = self._dg_dx_bl(x_vec)

        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        for _, k, l in self._nonzero_christoffels_list_bl[0:4]:
            val1 = dgdx[l, 0, k] + dgdx[k, 0, l]
            val2 = dgdx[l, 3, k] + dgdx[k, 3, l]
            chl[0, k, l] = chl[0, l, k] = 0.5 * (
                g_contra[0, 0] * (val1) + g_contra[0, 3] * (val2)
            )
            chl[3, k, l] = chl[3, l, k] = 0.5 * (
                g_contra[3, 0] * (val1) + g_contra[3, 3] * (val2)
            )
        for i, k, l in self._nonzero_christoffels_list_bl[8:16]:
            chl[i, k, l] = 0.5 * (
                g_contra[i, i] * (dgdx[l, i, k] + dgdx[k, i, l] - dgdx[i, k, l])
            )
        for i, k, l in self._nonzero_christoffels_list_bl[16:20]:
            chl[i, k, l] = chl[i, l, k] = 0.5 * (
                g_contra[i, i] * (dgdx[l, i, k] + dgdx[k, i, l] - dgdx[i, k, l])
            )

        return chl

    def _f_vec(self, lambda_, vec):
        """
        Returns f_vec for Kerr-Newman Metric in chosen coordinates
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
        if self.coords.system == "BoyerLindquist":
            return self._f_vec_bl(lambda_, vec)

        raise CoordinateError(
            "'f_vec' for Kerr-Newman Metric is available only in Boyer-Lindquist Coordinates."
        )

    def _f_vec_bl(self, lambda_, vec):
        """
        Returns f_vec for Kerr-Newman Metric \
        in Boyer-Lindquist Coordinates
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
            f_vec for Kerr-Newman Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)

        """
        chl = self.christoffels(vec[:4])
        F_contra = self.em_tensor_contravariant(vec[:4])
        g_cov = self.metric_covariant(vec[:4])

        vals = np.zeros(shape=vec.shape, dtype=vec.dtype)

        vals[:4] = vec[4:]

        vals[4] = -2.0 * (
            chl[0, 0, 1] * vec[4] * vec[5]
            + chl[0, 0, 2] * vec[4] * vec[6]
            + chl[0, 1, 3] * vec[5] * vec[7]
            + chl[0, 2, 3] * vec[6] * vec[7]
        )
        vals[5] = -1.0 * (
            chl[1, 0, 0] * vec[4] * vec[4]
            + 2 * chl[1, 0, 3] * vec[4] * vec[7]
            + chl[1, 1, 1] * vec[5] * vec[5]
            + 2 * chl[1, 1, 2] * vec[5] * vec[6]
            + chl[1, 2, 2] * vec[6] * vec[6]
            + chl[1, 3, 3] * vec[7] * vec[7]
        )
        vals[6] = -1.0 * (
            chl[2, 0, 0] * vec[4] * vec[4]
            + 2 * chl[2, 0, 3] * vec[4] * vec[7]
            + chl[2, 1, 1] * vec[5] * vec[5]
            + 2 * chl[2, 1, 2] * vec[5] * vec[6]
            + chl[2, 2, 2] * vec[6] * vec[6]
            + chl[2, 3, 3] * vec[7] * vec[7]
        )
        vals[7] = -2.0 * (
            chl[3, 0, 1] * vec[4] * vec[5]
            + chl[3, 0, 2] * vec[4] * vec[6]
            + chl[3, 1, 3] * vec[5] * vec[7]
            + chl[3, 2, 3] * vec[6] * vec[7]
        )

        vals[4:] -= self.q.value * (F_contra @ vec[4:] @ g_cov)

        return vals

    def em_potential_covariant(self, x_vec):
        """
        Returns Covariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        _, r, th, _ = x_vec
        M, a, Q = self.M.value, self.a.value, self.Q.value

        alpha = super().alpha(M, a)
        # Geometrized Charge
        r_Q = np.sqrt((Q**2 * _G * _Cc) / _c**4)
        rho2 = super().rho(r, th, M, a) ** 2

        A = np.zeros((4,), dtype=float)
        A[0] = r * r_Q / rho2
        A[3] = -r * alpha * r_Q * np.sin(th) ** 2 / rho2

        return A

    def em_potential_contravariant(self, x_vec):
        """
        Returns Contravariant Electromagnetic 4-Potential
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Contravariant Electromagnetic 4-Potential
            Numpy array of shape (4,)

        """
        A_cov = self.em_potential_covariant(x_vec)
        g_contra = self.metric_contravariant(x_vec)

        return g_contra @ A_cov

    def em_tensor_covariant(self, x_vec):
        """
        Returns Covariant Electromagnetic Tensor
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Electromagnetic Tensor
            Numpy array of shape (4, 4)

        """
        _, r, th, _ = x_vec
        M, a, Q = self.M.value, self.a.value, self.Q.value

        alpha = super().alpha(M, a)
        r_Q = np.sqrt((Q**2 * _G * _Cc) / _c**4)
        rho2 = super().rho(r, th, M, a) ** 2
        # Partial derivatives of rho2
        drho2_dr = 2 * r
        drho2_dtheta = -(alpha**2 * np.sin(2 * th))

        F = np.zeros((4, 4), dtype=float)

        F[0, 1] = -(r_Q * (rho2 - drho2_dr * r)) / (rho2**2)
        F[1, 0] = -F[0, 1]
        F[0, 2] = (r * r_Q * drho2_dtheta) / (rho2**2)
        F[2, 0] = -F[0, 2]
        F[1, 3] = (
            (1 / rho2**2) * (alpha * r_Q * np.sin(th) ** 2) * (rho2 - 2 * r**2)
        )
        F[3, 1] = -F[1, 3]
        F[2, 3] = (
            (1 / rho2**2)
            * (alpha * r_Q * r * np.sin(2 * th))
            * (rho2 + (alpha * np.sin(th)) ** 2)
        )
        F[3, 2] = -F[2, 3]

        return F

    def em_tensor_contravariant(self, x_vec):
        """
        Returns Contravariant Electromagnetic Tensor
        Specific to Kerr-Newman Geometries

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Contravariant Electromagnetic Tensor
            Numpy array of shape (4, 4)

        """
        F_cov = self.em_tensor_covariant(x_vec)
        g_contra = self.metric_contravariant(x_vec)

        F_contra = g_contra @ F_cov @ g_contra

        return F_contra
