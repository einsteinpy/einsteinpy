import numpy as np
from astropy import units as u

from einsteinpy import constant
from einsteinpy.metric import BaseMetric
from einsteinpy.utils import CoordinateError

_c = constant.c.value


class Kerr(BaseMetric):
    """
    Class for defining the Kerr Geometry

    """

    @u.quantity_input(M=u.kg, a=u.one)
    def __init__(self, coords, M, a):
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

        """
        super().__init__(
            coords=coords,
            M=M,
            a=a,
            name="Kerr Metric",
            metric_cov=self.metric_covariant,
            christoffels=self._christoffels,
            f_vec=self._f_vec,
        )
        # Precomputed list of tuples, containing indices \
        # of non-zero Christoffel Symbols for Kerr Metric \
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
        Returns Covariant Kerr Metric Tensor \
        in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr Metric Tensor in chosen Coordinates
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
            "Kerr Metric is available only in Boyer-Lindquist Coordinates."
        )

    def _g_cov_bl(self, x_vec):
        """
        Returns Covariant Kerr Metric Tensor \
        in Boyer-Lindquist coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Kerr Metric Tensor \
            in Boyer-Lindquist coordinates
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s, M, a = self.sch_rad, self.M.value, self.a.value
        alpha = super().alpha(M, a)
        sg, dl = super().sigma(r, th, M, a), super().delta(r, M, a)

        g_cov_bl = np.zeros(shape=(4, 4), dtype=float)

        g_cov_bl[0, 0] = (1 - (r_s * r / sg)) * _c ** 2
        g_cov_bl[1, 1] = -(sg / dl)
        g_cov_bl[2, 2] = -sg
        g_cov_bl[3, 3] = -(
            ((r ** 2) + (alpha ** 2) + ((r_s * r * (alpha * np.sin(th)) ** 2) / sg))
            * (np.sin(th) ** 2)
        )
        g_cov_bl[0, 3] = g_cov_bl[3, 0] = (_c * r_s * r * alpha * (np.sin(th) ** 2)) / (
            sg
        )

        return g_cov_bl

    def _dg_dx_bl(self, x_vec):
        """
        Returns derivative of each Kerr Metric component \
        w.r.t. coordinates in Boyer-Lindquist Coordinate System

        Parameters
        ----------
        x_vec : array_like
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
        r_s, M, a = self.sch_rad, self.M.value, self.a.value
        alpha = super().alpha(M, a)
        sg, dl = super().sigma(r, th, M, a), super().delta(r, M, a)

        dgdx = np.zeros(shape=(4, 4, 4), dtype=float)

        # Metric is invariant on t & phi
        # Differentiation of metric wrt r
        def due_to_r():
            nonlocal dgdx
            dsdr = 2 * r
            dddr = 2 * r - r_s
            tmp = r_s * (sg - r * dsdr) / (sg ** 2)  # r_s * d (r/sg) / dr
            dgdx[1, 0, 0] = -tmp * _c ** 2
            dgdx[1, 1, 1] = -(dsdr - (sg * (dddr / dl))) / dl
            dgdx[1, 2, 2] = -dsdr
            dgdx[1, 3, 3] = (-2 * r + ((alpha * np.sin(th)) ** 2) * tmp) * (
                np.sin(th) ** 2
            )
            dgdx[1, 0, 3] = dgdx[1, 3, 0] = _c * alpha * (np.sin(th) ** 2) * tmp

        # Differentiation of metric wrt theta
        def due_to_theta():
            nonlocal dgdx
            dsdth = -(alpha ** 2) * np.sin(2 * th)
            tmp = (
                ((_c / sg) ** 2) * r_s * r * dsdth
            )  # (- _c**2 * r_s * r) * d (1/sg) / dth
            dgdx[2, 0, 0] = tmp
            dgdx[2, 1, 1] = -(dsdth / dl)
            dgdx[2, 2, 2] = -dsdth
            dgdx[2, 3, 3] = -np.sin(2 * th) * ((r ** 2) + (alpha ** 2)) - (
                r_s * r * alpha ** 2
            ) * ((np.sin(th) / sg) ** 2) * (
                2 * sg * np.sin(2 * th) - (np.sin(th) ** 2) * dsdth
            )
            dgdx[2, 0, 3] = dgdx[2, 3, 0] = ((_c * alpha * r_s * r) / (sg ** 2)) * (
                sg * np.sin(2 * th) - dsdth * np.sin(th) ** 2
            )

        due_to_r()
        due_to_theta()

        return dgdx

    def _christoffels(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr Metric in chosen Coordinates

        Parameters
        ----------
        x_vec : array_like
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Kerr Metric \
            in chosen Coordinates
            Numpy array of shape (4,4,4)

        Raises
        ------
        CoordinateError
            Raised, if the Christoffel symbols are not \
            available in the supplied Coordinate System

        """
        if self.coords.system == "BoyerLindquist":
            return self._ch_sym_bl(x_vec)

        raise CoordinateError(
            "Christoffel Symbols for Kerr Metric are available only in Boyer-Lindquist Coordinates."
        )

    def _ch_sym_bl(self, x_vec):
        """
        Returns Christoffel Symbols for Kerr Metric \
        in Boyer-Lindquist Coordinates

        Parameters
        ----------
        x_vec : array_like
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
        Returns f_vec for Kerr Metric in chosen coordinates
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
            f_vec for Kerr Metric in chosen coordinates
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
            "'f_vec' for Kerr Metric is available only in Boyer-Lindquist Coordinates."
        )

    def _f_vec_bl(self, lambda_, vec):
        """
        Returns f_vec for Kerr Metric \
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
