# Issues:
# 1. Sign of f_vec
# 2. NullBundle


import warnings

import numpy as np

from einsteinpy.geodesic import Geodesic
from einsteinpy.geodesic.utils import v1  # - ????
from einsteinpy.integrators import RK45


class NullGeodesic(Geodesic):
    """
    Class for defining Null-like Geodesics

    """

    def __init__(
        self, metric, coords, end_lambda, step_size=1e-3, return_cartesian=True
    ):
        """
        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric, in which Geodesics are to be calculated
        coords : ~einsteinpy.coordinates.core.*
            Coordinate system, in which Metric is to be represented
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
        step_size : float, optional
            Size of each geodesic integration step
            Defaults to ``1e-3``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        """
        super().__init__(
            time_like=False,
            metric=metric,
            coords=coords,
            end_lambda=end_lambda,
            step_size=step_size,
            return_cartesian=return_cartesian,
        )

        self._state = self._calculate_state()
        # Trajectory integration will be done
        # by calculate_trajectory() in Geodesic

    def _calculate_state(self):
        """
        Prepares and returns the Initial State Vector of the massless test particle

        # Returns
        # -------
        # state : ~numpy.ndarray
        #     Initial State Vector of the massless test particle
        #     Length-8 Array

        # Raises
        # ------
        # TypeError
        #     If there is a mismatch between the coordinates class of ``self.coords`` and \
        #     coordinate class, with which ``self.metric`` object has been instantiated

        """
        if self.coords.system != self.metric.coords.system:
            raise TypeError(
                "Coordinate System Mismatch between Metric object and supplied initial coordinates."
            )

        x4 = self.coords.position()
        al = x4[1]
        be = self.metric.sch_rad * 5
        i = x4[2]

        g_cov_mat = self.metric.metric_covariant(x4)

        # CHECK UNITS - ?????
        E = 1.0

        L = -al * E * np.sqrt(1 - np.cos(i) ** 2)
        Q = (E ** 2) * (be ** 2 + (np.cos(i) ** 2) * (al ** 2 - 1))

        k_t = -E
        k_p = L
        k_th = np.sign(be) * np.sqrt(
            np.abs(Q - (L * np.cot(i)) ** 2 + (E * np.cos(i)) ** 2)
        )
        k_r = v1(g_cov_mat, k_t, k_th, k_p)

        v4 = [k_t, k_p, k_th, k_r]

        state = np.hstack((x4, v4))

        return state


class NullBundle:
    """
    Class for generating a photon sheet and performing
    geodesic integration for Radiative Transfer applications

    """
