import numpy as np

from einsteinpy import constant
from einsteinpy.metric import BaseMetric

_c = constant.c.value


class Kerr(BaseMetric):
    """
    Class for defining the Kerr Geometry

    """

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
            Spin Parameter

        """
        # Precomputed list of tuples, containing indices \
        # of non-zero Christoffel Symbols for Kerr Metric \
        # in Boyer-Lindquist Coordinates
        self.nonzero_christoffels_list_bl = [
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
            name="Kerr Metric",
            metric_cov=self.metric_covariant,
            christoffels=self._christoffels,
            f_vec=self._f_vec,
        )

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr Metric Tensor in chosen Coordinates
            Numpy array of shape (4,4)

        """
        if self.coords == "KS":
            return self._g_cov_ks(x_vec)

        return self._g_cov_bl(x_vec)

    def _g_cov_bl(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s, M, a, c2 = self.sch_rad, self.M, self.a, _c ** 2
        alpha = super().alpha(M, a)
        sg, dl = super().sigma(r, th, M, a), super().delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = 1 - (r_s * r / sg)
        g_cov_bl[1, 1] = (sg / dl) * (-1 / c2)
        g_cov_bl[2, 2] = -1 * sg / c2
        g_cov_bl[3, 3] = (
            (-1 / c2)
            * (
                (r ** 2)
                + (alpha ** 2)
                + (r_s * r * (np.sin(th) ** 2) * ((alpha ** 2) / sg))
            )
            * (np.sin(th) ** 2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = (
            r_s * r * alpha * (np.sin(th) ** 2) / (sg * _c)
        )

        return g_cov_bl

    def _g_cov_ks(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
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
        Returns derivative of each Kerr Metric component \
        w.r.t. coordinates in Boyer-Lindquist Coordinate System

        Parameters
        ----------
        x_vec : ~numpy.ndarray
                Position 4-Vector

        Returns
        -------
        dgdx : ~numpy.ndarray
            Array, containing derivative of each Kerr Metric \
            component w.r.t. coordinates \
            in Boyer-Lindquist Coordinate System
            Numpy array of shape (4,4,4)
            dgdx[0], dgdx[1], dgdx[2] & dgdx[3] contain \
            derivatives of metric w.r.t. t, r, theta & phi respectively

        """
        r, th = x_vec[1], x_vec[2]
        r_s, M, a, c2 = self.sch_rad, self.M, self.a, _c ** 2
        alpha = super().alpha(M, a)
        sg, dl = super().sigma(r, th, M, a), super().delta(r, M, a)

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
                * (2 * r + (alpha ** 2) * (np.sin(th) ** 2) * tmp)
                * (np.sin(th) ** 2)
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = (1 / _c) * (alpha * (np.sin(th) ** 2) * tmp)

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            dsdth = -2 * (alpha ** 2) * np.cos(th) * np.sin(th)
            tmp = (-1 / sg) * r_s * r * dsdth / sg
            dgdx[2, 0, 0] = -1 * tmp
            dgdx[2, 1, 1] = (-1 / c2) * (dsdth / dl)
            dgdx[2, 2, 2] = (-1 / c2) * dsdth
            dgdx[2, 3, 3] = (-1 / c2) * (
                2 * np.sin(th) * np.cos(th) * ((r ** 2) + (alpha ** 2))
                + tmp * (alpha ** 2) * (np.sin(th) ** 4)
                + (4 * (np.sin(th) ** 3) * np.cos(th) * (alpha ** 2) * r * r_s / sg)
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = (alpha / _c) * (
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
        Returns Christoffel Symbols for Kerr Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr Metric \
            in chosen Coordinates
            Numpy array of shape (4,4,4)

        """
        if self.coords == "KS":
            return self._ch_sym_ks(x_vec)

        return self._ch_sym_bl(x_vec)

    def _ch_sym_bl(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr Metric \
            in Boyer-Lindquist Coordinates
            Numpy array of shape (4,4,4)

        """
        g_contra = self.metric_contravariant(x_vec)
        dgdx = self._dg_dx_bl(x_vec)

        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        for _, k, l in self.nonzero_christoffels_list_bl[0:4]:
            val1 = dgdx[l, 0, k] + dgdx[k, 0, l]
            val2 = dgdx[l, 3, k] + dgdx[k, 3, l]
            chl[0, k, l] = chl[0, l, k] = 0.5 * (
                g_contra[0, 0] * (val1) + g_contra[0, 3] * (val2)
            )
            chl[3, k, l] = chl[3, l, k] = 0.5 * (
                g_contra[3, 0] * (val1) + g_contra[3, 3] * (val2)
            )
        for i, k, l in self.nonzero_christoffels_list_bl[8:16]:
            chl[i, k, l] = 0.5 * (
                g_contra[i, i] * (dgdx[l, i, k] + dgdx[k, i, l] - dgdx[i, k, l])
            )
        for i, k, l in self.nonzero_christoffels_list_bl[16:20]:
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
        Returns f_vec for Kerr Metric in chosen coordinates
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
            f_vec for Kerr Metric in chosen coordinates
            Numpy array of shape (8)

        """
        if self.coords == "KS":
            return self._f_vec_ks(lambda_, x_vec)

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

        vec : ~numpy.ndarray
            Length-8 Vector, containing 4-Position & 4-Velocity

        Returns
        -------
        ~numpy.ndarray
            f_vec for Kerr Metric in Boyer-Lindquist Coordinates
            Numpy array of shape (8)

        """
        chl = self.christoffels(vec[:4])
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

        vec : ~numpy.ndarray
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

            for index in range(4):
                chl[i, j, k] |= g_contra[i, index] & (
                    dgdx[k, index, j] | dgdx[j, index, k] | dgdx[index, j, k]
                )

            if chl[i, j, k]:
                vcl.append((i, j, k))

        return vcl

    # Hiding unrelated methods
    em_potential_covariant = BaseMetric._private
    em_potential_contravariant = BaseMetric._private
    em_tensor_covariant = BaseMetric._private
    em_tensor_contravariant = BaseMetric._private
