import numpy as np

from einsteinpy import constant
from einsteinpy.metric import BaseMetric


__all__ = ['Schwarzschild']


_c = constant.c.value


class Schwarzschild(BaseMetric):
    """
    Class for defining Schwarzschild Geometry

    """

    def __init__(self, M):
        """
        Constructor

        Parameters
        ----------
        M : float
            Mass of gravitating body, e.g. Black Hole

        """
        super().__init__(
            coords="S",
            M=M,
            name="Schwarzschild Metric",
            metric_cov=self.metric_covariant,
            christoffels=self._christoffels,
            f_vec=self._f_vec,
        )

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Schwarzschild Metric Tensor \
        in Schwarschild Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Covariant Schwarzschild Metric Tensor
            Numpy array of shape (4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s = self.sch_rad
        g_cov = np.zeros(shape=(4, 4), dtype=float)

        tmp, c2 = 1.0 - (r_s / r), _c ** 2
        g_cov[0, 0] = tmp
        g_cov[1, 1] = -1.0 / (tmp * c2)
        g_cov[2, 2] = -1 * (r ** 2) / c2
        g_cov[3, 3] = -1 * ((r * np.sin(th)) ** 2) / c2

        return g_cov

    def _christoffels(self, x_vec):
        """
        Returns the Christoffel Symbols for Schwarzschild Metric
        in Schwarzschild Coordinates

        Parameters
        ----------
        x_vec : ~numpy.ndarray
            Position 4-Vector

        Returns
        -------
        ~numpy.ndarray
            Christoffel Symbols for Schwarzschild Metric \
            in Schwarzschild Coordinates
            Numpy array of shape (4,4,4)

        """
        r, th = x_vec[1], x_vec[2]
        r_s = self.sch_rad
        c2 = _c ** 2
        chl = np.zeros(shape=(4, 4, 4), dtype=float)

        chl[1, 0, 0] = 0.5 * r_s * (r - r_s) * c2 / (r ** 3)
        chl[1, 1, 1] = 0.5 * r_s / (r_s * r - r ** 2)
        chl[1, 2, 2] = r_s - r
        chl[1, 3, 3] = (r_s - r) * (np.sin(th) ** 2)
        chl[0, 0, 1] = chl[0, 1, 0] = -chl[1, 1, 1]
        chl[2, 2, 1] = chl[2, 1, 2] = chl[3, 3, 1] = chl[3, 1, 3] = 1 / r
        chl[2, 3, 3] = -np.cos(th) * np.sin(th)
        chl[3, 3, 2] = chl[3, 2, 3] = 1 / np.tan(th)

        return chl

    def _f_vec(self, lambda_, vec):
        """
        Returns f_vec for Schwarzschild Metric
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
            f_vec for Schwarzschild Metric
            Numpy array of shape (8)

        """
        chl = self.christoffels(vec[:4])
        vals = np.zeros(shape=vec.shape, dtype=vec.dtype)

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
