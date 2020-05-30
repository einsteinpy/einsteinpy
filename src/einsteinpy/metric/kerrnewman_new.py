import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.metric import Metric


_c = constant.c.value

# Precomputed list of tuples, containing indices \
# of non-zero Christoffel Symbols for Kerr-Newman Metric
nonzero_christoffels_list = [
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

class KerrNewman(Metric):
    """
    Class for defining Kerr-Newman Goemetry
    """

    @u.quantity_input(M=u.kg, a=u.m, Q=u.C)
    def __init__(self, coords, M, a, Q):
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
            Spin Parameter, 0 <= a <= 1 - ???? (Only if Geom Units)

        """
        super().__init__(
            coords = coords,
            M = M,
            a = a,
            Q = Q,
            name = "Kerr-Newman Metric",
            metric_cov = self.metric_covariant,
            christoffels = self.christoffels_,
            f_vec = self.f_vec_
        )

    # - ????? (Overrides Metric.metric_covariant() - ✅)
    def metric_covariant(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Kerr-Newman Metric Tensor in chosen Coordinates
            Numpy array of shape (4,4)

        """
        if self.coords == "BL":
            return self._g_cov_bl(x_vec)
        
        elif self.coords == "KS":
            return self._g_cov_ks(x_vec)
        
        # Default choice
        return self._g_cov_bl(x_vec)
    
    # - ????? (Contravariant form returned by super class) - ✅

    def _g_cov_bl(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Kerr-Newman Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        M, a, c2 = self.M, self.a, _c ** 2
        rho2, dl = Metric.rho(r, th, a), Metric.delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = (dl - ((a * np.sin(th)) ** 2)) / (rho2)
        g_cov_bl[1, 1] = -rho2 / (dl * c2)
        g_cov_bl[2, 2] = -rho2 / c2
        g_cov_bl[3, 3] = (
            (((a * np.sin(th)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2))
            * (np.sin(th) ** 2)
            / (rho2 * c2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = (
            -a * (np.sin(th) ** 2) * (dl - (r ** 2) - (a ** 2)) / (rho2 * _c)
        )

        return g_cov_bl

    def _g_cov_ks(self, x_vec):
        """
        Returns Covariant Kerr-Newman Metric Tensor \
        in Kerr-Schild coordinates

        Parameters
        ----------
        x_vec : numpy.array
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
        x_vec : numpy.array
                Position 4-Vector
        
        Returns
        -------
        dgdx : ~numpy.array
            Array, containing derivative of each Kerr-Newman \
            Metric component w.r.t. coordinates \
            in Boyer-Lindquist Coordinate System
            Numpy array of shape (4,4,4)
            dgdx[0], dgdx[1], dgdx[2] & dgdx[3] contain \
            derivatives of metric w.r.t. t, r, theta & phi respectively

        """
        r, th = x_vec[1], x_vec[2]
        M, a, c2 = self.M, self.a, _c ** 2
        rho2, dl = Metric.rho(r, th, a), Metric.delta(r, M, a)

        dgdx = np.zeros(shape=(4, 4, 4), dtype=float)
        
        # Metric is invariant on t & phi
        # Differentiation of metric wrt r
        def due_to_r():
            nonlocal dgdx
            drh2dr = 2 * r
            dddr = 2 * r - self.sch_rad.value
            dgdx[1, 0, 0] = (dddr * rho2 - drh2dr * (dl - (a * np.sin(th)) ** 2)) / (
                rho2 ** 2
            )
            dgdx[1, 1, 1] = (-1 / (c2 * (dl ** 2))) * (drh2dr * dl - dddr * rho2)
            dgdx[1, 2, 2] = -drh2dr / c2
            dgdx[1, 3, 3] = ((np.sin(th) ** 2) / (c2 * (rho2 ** 2))) * (
                (
                    (((a * np.sin(th)) ** 2) * dddr - 4 * (r ** 3) - 4 * (r * (a ** 2)))
                    * rho2
                )
                - (drh2dr * (((a * np.sin(th)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2)))
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = (
                (-a) * (np.sin(th) ** 2) / (_c * (rho2 ** 2))
            ) * ((dddr - 2 * r) * rho2 - drh2dr * (dl - r ** 2 - a ** 2))

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            drh2dth = -2 * (a ** 2) * np.cos(th) * np.sin(th)
            dgdx[2, 0, 0] = (
                (-2 * (a ** 2) * np.sin(th) * np.cos(th)) * rho2
                - drh2dth * (dl - ((a * np.sin(th)) ** 2))
            ) / (rho2 ** 2)
            dgdx[2, 1, 1] = -drh2dth / (c2 * dl)
            dgdx[2, 2, 2] = -drh2dth / c2
            dgdx[2, 3, 3] = (1 / (c2 * (rho2 ** 2))) * (
                (
                    (
                        (4 * (a ** 2) * (np.sin(th) ** 3) * np.cos(th) * dl)
                        - (2 * np.sin(th) * np.cos(th) * ((r ** 2 + a ** 2) ** 2))
                    )
                    * rho2
                )
                - (
                    drh2dth
                    * (((a * np.sin(th)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2))
                    * (np.sin(th) ** 2)
                )
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = (
                (-a * (dl - r ** 2 - a ** 2)) / (_c * (rho2 ** 2))
            ) * (
                (2 * np.sin(th) * np.cos(th) * rho2) - (drh2dth * (np.sin(th) ** 2))
            )

        due_to_r()
        due_to_theta()

        return dgdx

    def _dg_dx_ks(self, x_vec):
        """
        Returns derivative of each Kerr-Newman Metric component \
        w.r.t. coordinates in Kerr-Schild Coordinate System

        Parameters
        ----------
        x_vec : numpy.array
                Position 4-Vector
        
        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

    def christoffels_(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Christoffel Symbols for Kerr-Newman \
            Metric in chosen Coordinates
            Numpy array of shape (4,4,4)
        
        """
        if self.coords == "BL":
            return self._ch_sym_bl(x_vec)
        
        elif self.coords == "KS":
            return self._ch_sym_ks(x_vec)
        
        # Default choice
        return self._ch_sym_bl(x_vec)

    def _ch_sym_bl(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr-Newman Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Christoffel Symbols for Kerr-Newman Metric \
            in Boyer-Lindquist Coordinates
            Numpy array of shape (4,4,4)
        
        """
        g_contra = self.metric_contravariant(x_vec)
        dgdx = self._dg_dx_bl(x_vec)

        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        for _, k, l in nonzero_christoffels_list[0:4]:
            val1 = dgdx[l, 0, k] + dgdx[k, 0, l]
            val2 = dgdx[l, 3, k] + dgdx[k, 3, l]
            chl[0, k, l] = chl[0, l, k] = 0.5 * (g_contra[0, 0] * (val1) + g_contra[0, 3] * (val2))
            chl[3, k, l] = chl[3, l, k] = 0.5 * (g_contra[3, 0] * (val1) + g_contra[3, 3] * (val2))
        for i, k, l in nonzero_christoffels_list[8:16]:
            chl[i, k, l] = 0.5 * (
                g_contra[i, i] * (dgdx[l, i, k] + dgdx[k, i, l] - dgdx[i, k, l])
            )
        for i, k, l in nonzero_christoffels_list[16:20]:
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
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates

        """
        # To be implemented after KS Coordinates
        raise NotImplementedError

    def f_vec_(self, lambda_, x_vec):
        """
        Returns f_vec for Kerr-Newman Metric in chosen coordinates
        To be used in solving for Geodesics

        Parameters
        ----------
        lambda_ : float
            Parameterizes current integration step
            Used by ODE Solver

        vec : numpy.array
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.array
            f_vec for Kerr-Newman Metric in chosen coordinates
            Numpy array of shape (8)
        
        """
        if self.coords == "BL":
            return self._f_vec_bl(lambda_, x_vec)
        
        elif self.coords == "KS":
            return self._f_vec_ks(lambda_, x_vec)
        
        # Default choice
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

        vec : numpy.array
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.array
            f_vec for Kerr-Newman Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)
        
        """
        chl = self.christoffels(
            vec[1], vec[2], self.M.value, self.a.value, self.Q.value
        )
        F_contra = self.maxwell_tensor_contravariant(
            vec[1], vec[2], self.a.value, self.Q.value, self.M.value
        )
        x_vec = [0, vec[1], vec[2], 0] # t & phi have no bearing on Metric
        g_cov = self.metric_covariant(x_vec)

        vals = np.zeros(shape=(8,), dtype=float)

        for i in range(4):
            vals[i] = vec[i + 4]
        
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
        vals[4:] -= self.Q.value * np.dot(
            vec[4:].reshape((4,)), g_cov @ F_contra
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

        vec : numpy.array
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        NotImplementedError
            To be implemented after KS Coordinates
        
        """
        # To be implemented after KS Coordinates
        raise NotImplementedError