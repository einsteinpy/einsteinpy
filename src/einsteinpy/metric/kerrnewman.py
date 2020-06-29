import numpy as np

from einsteinpy import constant
from einsteinpy.metric import BaseMetric

_c = constant.c.value
_G = constant.G.value
_Cc = constant.coulombs_const.value


class KerrNewman(BaseMetric):
    """
    Class for defining Kerr-Newman Goemetry

    """

    def __init__(self, coords, M, a, Q, q=0.0):
        """
        Constructor

        Parameters
        ----------
        coords : string
            Coordinate system, in which Metric is to be represented
            "BL" - Boyer-Lindquist: Applicable to Kerr-Newman solutions
            "KS" - Kerr-Schild: Useful for adding perturbations to Kerr-Newman solutions
        M : float
            Mass of gravitating body, e.g. Black Hole
        a : float
            Spin Parameter
        Q : float
            Charge on gravitating body, e.g. Black Hole
        q : float, optional
            Charge, per unit mass, of the test particle
            Defaults to ``O``

        """
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

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr-Newman Metric Tensor in chosen Coordinates
            Numpy array of shape (4,4)

        """
        if self.coords == "KS":
            return self._g_cov_ks(x_vec)

        return self._g_cov_bl(x_vec)

    def _g_cov_bl(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr-Newman Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        M, a = self.M, self.a
        alpha = super().alpha(M, a)
        rho2, del_ = super().rho(r, th, M, a) ** 2, super().delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = (_c ** 2) * ((del_ - ((alpha * np.sin(th)) ** 2)) / (rho2))
        g_cov_bl[1, 1] = -rho2 / del_
        g_cov_bl[2, 2] = -rho2
        g_cov_bl[3, 3] = -(
            (np.sin(th) ** 2)
            * (((r ** 2 + alpha ** 2) ** 2 - del_ * (alpha * np.sin(th)) ** 2) / rho2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = _c * (
            (-alpha * (np.sin(th) ** 2) * (del_ - (r ** 2) - (alpha ** 2))) / rho2
        )

        return g_cov_bl

    def _g_cov_ks(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in Kerr-Schild coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

    def _dg_dx_bl(self, x_vec):
        """
        Returns derivative of each Kerr-Newman Metric component \
        w.r.t. coordinates in Boyer-Lindquist Coordinate System

        Parameters
        ----------
        x_vec : ~numpy.ndarray
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
        M, a, c2 = self.M, self.a, _c ** 2
        alpha = super().alpha(M, a)
        rho2, del_ = super().rho(r, th, M, a) ** 2, super().delta(r, M, a)

        dgdx = np.zeros(shape=(4, 4, 4), dtype=float)

        # Metric is invariant on t & phi
        # Differentiation of metric wrt r
        def due_to_r():
            nonlocal dgdx
            drh2dr = 2 * r
            dddr = 2 * r - self.sch_rad
            dgdx[1, 0, 0] = (
                dddr * rho2 - drh2dr * (del_ - (alpha * np.sin(th)) ** 2)
            ) / (rho2 ** 2)
            dgdx[1, 1, 1] = (-1 / (c2 * (del_ ** 2))) * (drh2dr * del_ - dddr * rho2)
            dgdx[1, 2, 2] = -drh2dr / c2
            dgdx[1, 3, 3] = ((np.sin(th) ** 2) / (c2 * (rho2 ** 2))) * (
                (
                    (
                        ((alpha * np.sin(th)) ** 2) * dddr
                        - 4 * (r ** 3)
                        - 4 * (r * (alpha ** 2))
                    )
                    * rho2
                )
                - (
                    drh2dr
                    * (
                        ((alpha * np.sin(th)) ** 2) * del_
                        - ((r ** 2 + alpha ** 2) ** 2)
                    )
                )
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = (
                (-alpha) * (np.sin(th) ** 2) / (_c * (rho2 ** 2))
            ) * ((dddr - 2 * r) * rho2 - drh2dr * (del_ - r ** 2 - alpha ** 2))

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            drh2dth = -2 * (alpha ** 2) * np.cos(th) * np.sin(th)
            dgdx[2, 0, 0] = (
                (-2 * (alpha ** 2) * np.sin(th) * np.cos(th)) * rho2
                - drh2dth * (del_ - ((alpha * np.sin(th)) ** 2))
            ) / (rho2 ** 2)
            dgdx[2, 1, 1] = -drh2dth / (c2 * del_)
            dgdx[2, 2, 2] = -drh2dth / c2
            dgdx[2, 3, 3] = (1 / (c2 * (rho2 ** 2))) * (
                (
                    (
                        (4 * (alpha ** 2) * (np.sin(th) ** 3) * np.cos(th) * del_)
                        - (2 * np.sin(th) * np.cos(th) * ((r ** 2 + alpha ** 2) ** 2))
                    )
                    * rho2
                )
                - (
                    drh2dth
                    * (
                        ((alpha * np.sin(th)) ** 2) * del_
                        - ((r ** 2 + alpha ** 2) ** 2)
                    )
                    * (np.sin(th) ** 2)
                )
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = (
                (-alpha * (del_ - r ** 2 - alpha ** 2)) / (_c * (rho2 ** 2))
            ) * ((2 * np.sin(th) * np.cos(th) * rho2) - (drh2dth * (np.sin(th) ** 2)))

        due_to_r()
        due_to_theta()

        return dgdx

    def _dg_dx_ks(self, x_vec):
        """
        Returns derivative of each Kerr-Newman Metric component \
        w.r.t. coordinates in Kerr-Schild Coordinate System

        Parameters
        ----------
        x_vec : ~numpy.ndarray
                Position 4-Vector

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

    def _christoffels(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr-Newman \
            Metric in chosen Coordinates
            Numpy array of shape (4,4,4)

        """
        if self.coords == "KS":
            return self._ch_sym_ks(x_vec)

        return self._ch_sym_bl(x_vec)

    def _ch_sym_bl(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
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

    def _ch_sym_ks(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric \
        in Kerr-Schild Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

    def _f_vec(self, lambda_, x_vec):
        """
        Returns f_vec for Kerr-Newman Metric in chosen coordinates
        To be used in solving for Geodesics

        Parameters
        ----------
        lambda_ : float
            Parameterizes current integration step
            Used by ODE Solver

        vec : ~numpy.ndarray
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.ndarray
            f_vec for Kerr-Newman Metric in chosen coordinates
            Numpy array of shape (8)

        """
        if self.coords == "KS":
            return self._f_vec_ks(lambda_, x_vec)

        return self._f_vec_bl(lambda_, x_vec)

    def _f_vec_bl(self, lambda_, vec):
        """
        Returns f_vec for Kerr-Newman Metric \
        in Boyer-Lindquist Coordinates
        To be used in solving for Geodesics

        Parameters
        ----------
        lambda_ : float
            Parameterizes current integration step
            Used by ODE Solver

        vec : ~numpy.ndarray
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.ndarray
            f_vec for Kerr-Newman Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)

        """
        chl = self.christoffels(vec[:4])
        F_contra = self.em_tensor_contravariant(vec[1], vec[2], self.M, self.a, self.Q)
        x_vec = np.array(
            [0, vec[1], vec[2], 0], dtype=float
        )  # t & phi have no bearing on Metric
        g_cov = self.metric_covariant(x_vec)

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

        vals[4:] -= (self.q * np.dot(vec[4:].reshape((4,)), g_cov @ F_contra)).reshape(
            4, 1
        )

        return vals

    def _f_vec_ks(self, lambda_, vec):
        """
        Returns f_vec for Kerr-Newman Metric \
        in Kerr-Schild Coordinates
        To be used in solving for Geodesics

        Parameters
        ----------
        lambda_ : float
            Parameterizes current integration step
            Used by ODE Solver

        vec : ~numpy.ndarray
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

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
        alpha = super().alpha(M, a)
        # Geometrized Charge
        r_Q = np.sqrt((Q ** 2 * _G * _Cc) / _c ** 4)
        rho2 = super().rho(r, theta, M, a) ** 2

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
        alpha = super().alpha(M, a)
        r_Q = np.sqrt((Q ** 2 * _G * _Cc) / _c ** 4)
        rho2 = super().rho(r, theta, M, a) ** 2
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
