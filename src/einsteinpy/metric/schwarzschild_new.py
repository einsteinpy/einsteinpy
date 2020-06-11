import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.metric import Metric


_c = constant.c.value

class Schwarzschild(Metric):
    """
    Class for defining Schwarzschild Geometry

    """

    @u.quantity_input(M=u.kg)
    def __init__(self, M):
        """
        Constructor

        Parameters
        ----------
        M : float
            Mass of gravitating body, e.g. Black Hole

        """
        super().__init__(
            coords = "S",
            M = M,
            name = "Schwarzschild Metric",
            metric_cov = self.metric_covariant,
            christoffels = self.christoffels_,
            f_vec = self.f_vec_
        )

    def metric_covariant(self, x_vec):
        """
        Returns Covariant Schwarzschild Metric Tensor \
        in Schwarschild Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Covariant Schwarzschild Metric Tensor
            Numpy array of shape (4,4)

        """
        # Since, we only have Schwarzschild coordinates to define this Tensor
        # this `metric_covariant()` function does not switch coordinates,
        # like in `Kerr` or `KerrNewman`
        r, th = x_vec[1], x_vec[2]
        r_s = self.sch_rad.value
        g_cov = np.zeros(shape=(4, 4), dtype=float)

        tmp, c2 = 1. - (r_s / r), _c ** 2
        g_cov[0, 0] = tmp
        g_cov[1, 1] = -1.0 / (tmp * c2)
        g_cov[2, 2] = -1 * (r ** 2) / c2
        g_cov[3, 3] = -1 * ((r * np.sin(th)) ** 2) / c2

        return g_cov

    # - ????? (Contravariant form returned by super class) - ✅

    def christoffels_(self, x_vec):
        """
        Returns the Christoffel Symbols for Schwarzschild Metric
        in Schwarzschild Coordinates

        Parameters
        ----------
        x_vec : numpy.array
            Position 4-Vector

        Returns
        -------
        ~numpy.array
            Christoffel Symbols for Schwarzschild Metric \
            in Schwarzschild Coordinates
            Numpy array of shape (4,4,4)

        """
        # Removing c, G as parameters - ?????
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

    def f_vec_(self, lambda_, vec):
        """
        Returns f_vec for Schwarzschild Metric
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
            f_vec for Schwarzschild Metric
            Numpy array of shape (8)
        
        """
        chl = self.christoffels(vec)
        vals = np.zeros(shape=(8,), dtype=float)

        vals[:4] = vec[4:8]
        vals[4] = -2 * chl[0, 0, 1] * vec[4] * vec[5]
        vals[5] = -1 * (
            chl[1, 0, 0] * (vec[4] ** 2)
            + chl[1, 1, 1] * (vec[5] ** 2)
            + chl[1, 2, 2] * (vec[6] ** 2)
            + chl[1, 3, 3] * (vec[7] ** 2)
        )
        vals[6] = -2 * chl[2, 2, 1] * vec[6] * vec[5] - 1 * \
            chl[2, 3, 3] * (vec[7] ** 2)
        
        vals[7] = -2 * ( chl[3, 3, 1] * vec[7] * vec[5] + 
            chl[3, 3, 2] * vec[7] * vec[6] 
        )

        return vals

# time_velocity to be moved to `coordinates`
# calculate_trajectory to be moved to `geodesic`