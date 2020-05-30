import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.metric import Metric

_c = constant.c.value

# Precomputed list of tuples, containing indices \
# of non-zero Christoffel Symbols for Kerr Metric
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


class Kerr(Metric):
    """
    Class for defining the Kerr Geometry
    """

    @u.quantity_input(M=u.kg, a=u.km)
    def __init__(self, coords, M, a):
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
            coords=coords,
            M=M,
            a=a,
            name="Kerr Metric",
            metric_cov=self.metric_covariant,
            christoffels=self.christoffels_,
            f_vec=self.f_vec_,
        )

    # - ????? (Overrides Metric.metric_covariant() - ✅)
    def metric_covariant(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Kerr Metric Tensor in chosen Coordinates
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
        Returns Covariant Kerr Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Kerr Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s, M, a, c2 = self.sch_rad.value, self.M, self.a, _c ** 2
        sg, dl = Metric.sigma(r, th, a), Metric.delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = 1 - (r_s * r / sg)
        g_cov_bl[1, 1] = (sg / dl) * (-1 / c2)
        g_cov_bl[2, 2] = -1 * sg / c2
        g_cov_bl[3, 3] = (
            (-1 / c2)
            * ((r ** 2) + (a ** 2) + (r_s * r * (np.sin(th) ** 2) * ((a ** 2) / sg)))
            * (np.sin(th) ** 2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = r_s * r * a * (np.sin(th) ** 2) / (sg * _c)

        return g_cov_bl

    def _g_cov_ks(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
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
        Returns derivative of each Kerr Metric component \
        w.r.t. coordinates in Boyer-Lindquist Coordinate System

        Parameters
        ----------
        x_vec : numpy.array
                Position 4-Vector
        
        Returns
        -------
        dgdx : ~numpy.array
            Array, containing derivative of each Kerr Metric \
            component w.r.t. coordinates \
            in Boyer-Lindquist Coordinate System
            Numpy array of shape (4,4,4)
            dgdx[0], dgdx[1], dgdx[2] & dgdx[3] contain \
            derivatives of metric w.r.t. t, r, theta & phi respectively

        """
        r, th = x_vec[1], x_vec[2]
        r_s, M, a, c2 = self.sch_rad.value, self.M, self.a, _c ** 2
        sg, dl = Metric.sigma(r, th, a), Metric.delta(r, M, a)

        dgdx = np.zeros(shape=(4, 4, 4), dtype=float)

        # Metric is invariant on t & phi
        # Differentiation of metric wrt r
        def due_to_r():
            nonlocal dgdx
            dsdr = 2 * r
            dddr = 2 * r - r_s
            tmp = (r_s * (sg - r * dsdr) / sg) * (1 / sg)
            dgdx[1, 0, 0] = -1 * tmp
            dgdx[1, 1, 1] = (-1 / c2) * (dsdr - (sg * (dddr / dl))) / dl
            dgdx[1, 2, 2] = (-1 / c2) * dsdr
            dgdx[1, 3, 3] = (
                (-1 / c2)
                * (2 * r + (a ** 2) * (np.sin(th) ** 2) * tmp)
                * (np.sin(th) ** 2)
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = (1 / _c) * (a * (np.sin(th) ** 2) * tmp)

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            dsdth = -2 * (a ** 2) * np.cos(th) * np.sin(th)
            tmp = (-1 / sg) * r_s * r * dsdth / sg
            dgdx[2, 0, 0] = -1 * tmp
            dgdx[2, 1, 1] = (-1 / c2) * (dsdth / dl)
            dgdx[2, 2, 2] = (-1 / c2) * dsdth
            dgdx[2, 3, 3] = (-1 / c2) * (
                2 * np.sin(th) * np.cos(th) * ((r ** 2) + (a ** 2))
                + tmp * (a ** 2) * (np.sin(th) ** 4)
                + (4 * (np.sin(th) ** 3) * np.cos(th) * (a ** 2) * r * r_s / sg)
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = (a / _c) * (
                (np.sin(th) ** 2) * tmp + (2 * np.sin(th) * np.cos(th) * r_s * r / sg)
            )

        due_to_r()
        due_to_theta()

        return dgdx

    def _dg_dx_ks(self, x_vec):
        """
        Returns derivative of each Kerr Metric component \
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
        Returns Christoffel Symbols for Kerr Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Christoffel Symbols for Kerr Metric \
            in chosen Coordinates
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
        Returns Christoffel Symbols for Kerr Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Christoffel Symbols for Kerr Metric \
            in Boyer-Lindquist Coordinates
            Numpy array of shape (4,4,4)
        
        """
        g_contra = self.metric_contravariant(x_vec)
        dgdx = self._dg_dx_bl(x_vec)

        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        for _, k, l in nonzero_christoffels_list[0:4]:
            val1 = dgdx[l, 0, k] + dgdx[k, 0, l]
            val2 = dgdx[l, 3, k] + dgdx[k, 3, l]
            chl[0, k, l] = chl[0, l, k] = 0.5 * (
                g_contra[0, 0] * (val1) + g_contra[0, 3] * (val2)
            )
            chl[3, k, l] = chl[3, l, k] = 0.5 * (
                g_contra[3, 0] * (val1) + g_contra[3, 3] * (val2)
            )
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
        Returns Christoffel Symbols for Kerr Metric \
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
        Returns f_vec for Kerr Metric in chosen coordinates
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
            f_vec for Kerr Metric in chosen coordinates
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
        Returns f_vec for Kerr Metric \
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
            f_vec for Kerr Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)
        
        """
        chl = self.christoffels(vec)
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

        return vals

    def _f_vec_ks(self, lambda_, vec):
        """
        Returns f_vec for Kerr Metric \
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

    @staticmethod
    def nonzero_christoffels():
        """
        Returns a list of tuples consisting of indices \
        of non-zero Christoffel Symbols in Kerr Metric, \
        computed in real-time

        Returns
        -------
        list
            List of tuples
            Each tuple (i,j,k) represents Christoffel Symbols, \
            with i as upper index and j, k as lower indices.

        """
        # Below is the code for algorithmically calculating
        # the indices of nonzero christoffel symbols in Kerr Metric.
        g_contra = np.zeros(shape=(4, 4), dtype=bool)
        dgdx = np.zeros(shape=(4, 4, 4), dtype=bool)

        g_contra[3, 0] = g_contra[0, 3] = True
        dgdx[1, 3, 0] = dgdx[1, 0, 3] = True
        dgdx[2, 0, 3] = dgdx[2, 3, 0] = True
        # Code Climate cyclomatic complexity hack ; Instead of "for i in range(4)"
        g_contra[0, 0] = g_contra[1, 1] = g_contra[2, 2] = g_contra[3, 3] = True
        dgdx[1, 0, 0] = dgdx[1, 1, 1] = dgdx[1, 2, 2] = dgdx[1, 3, 3] = True
        dgdx[2, 0, 0] = dgdx[2, 1, 1] = dgdx[2, 2, 2] = dgdx[2, 3, 3] = True
        # hack ends

        chl = np.zeros(shape=(4, 4, 4), dtype=bool)
        tmp = np.array([i for i in range(4 ** 3)])
        vcl = list()

        for t in tmp:
            i = int(t / (4 ** 2)) % 4
            j = int(t / 4) % 4
            k = t % 4

            for l in range(4):
                chl[i, j, k] |= g_contra[i, l] & (
                    dgdx[k, l, j] | dgdx[j, l, k] | dgdx[l, j, k]
                )

            if chl[i, j, k]:
                vcl.append((i, j, k))

        return vcl
