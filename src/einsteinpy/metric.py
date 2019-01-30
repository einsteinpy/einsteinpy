import numpy as np

import astropy.units as u

from einsteinpy import constant
from einsteinpy.utils import schwarzschild_radius

class SchwarzschildMetric:

    def __init__(self, t, r, theta, phi, M):
        """
        Parameters
        ----------
        r : float

        """
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
        self.M = M

    @classmethod
    def from_position(cls, t, r, theta, phi, M):
        return cls(t, r, theta, phi, M)

    def _Rused(self):
        return schwarzschild_radius(self.M)/self.r

    @property
    def g00(self):
        return 1 - self._Rused()

    @property
    def g11(self):
        return -1 * (1/(1 - self._Rused()))

    @property
    def g22(self):
        return -1 * (self.r ** 2)

    @property
    def g33(self):
        return -1 * (self.r ** 2) * (np.sin(self.theta)**2)

    def print_metric(self):
        metric = np.zeros((4,4), dtype=float)
        metric[0][0] = self.g00
        metric[1][1] = self.g11
        metric[2][2] = self.g22
        metric[3][3] = self.g33
        return metric
