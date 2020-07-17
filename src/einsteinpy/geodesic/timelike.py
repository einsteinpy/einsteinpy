import warnings

import numpy as np

from einsteinpy.geodesic import Geodesic
from einsteinpy.integrators import RK45


class Timelike(Geodesic):
    """
    Class for defining Time-like Geodesics

    """

    def __init__(
        self, metric, coords, end_lambda, step_size=1e-3, return_cartesian=True
    ):
        """
        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric, in which Geodesics are to be calculated
        coords : ~einsteinpy.coordinates.differential.*
            Coordinate system, in which Metric is to be represented
        end_lambda : float
            Affine Parameter, Lambda, where iterations will stop
            Equivalent to Proper Time for Timelike Geodesics
        step_size : float, optional
            Size of each geodesic integration step
            Defaults to ``1e-3``
        return_cartesian : bool, optional
            Whether to return calculated values in Cartesian Coordinates
            Defaults to ``True``

        """
        super().__init__(
            time_like=True,
            metric=metric,
            coords=coords,
            end_lambda=end_lambda,
            step_size=step_size,
            return_cartesian=return_cartesian,
        )
