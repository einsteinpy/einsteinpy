import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.utils import scalar_factor, schwarzschild_radius


class RobertsonWalkerFlat:
    """
    Class for defining a Robertson Walker Flat Metric
    """

    def __init__(self, t, x, y, z, era="md"):
        """
        Constructor.

        Parameters
        ----------
        t : float
            Time for the event
        x : float
            Spatial Cartesian coordinates
        y : float
            Spatial Cartesian coordinates
        z : float
            Spatial Cartesian coordinates
        era : string, optional
            Can be chosen from md (Matter Dominant),
            rd (Radiation Dominant) and ded (Dark Energy Dominant)

        """
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self.era = era

    @classmethod
    @u.quantity_input(x=u.km, y=u.km, z=u.km)
    def from_position(cls, t, x, y, z, era):
        return cls(t, x, y, z, era)

    @property
    def g00(self):
        return -1

    @property
    def g11(self):
        return -1 * (1 / (1 - self._rused()))

    @property
    def g22(self):
        return -1 * (self.r ** 2)

    @property
    def g33(self):
        return -1 * (self.r ** 2) * (np.sin(self.theta) ** 2)

    @property
    def christ_sym1_00(self):
        num1 = (-2 * constant.G * self.M) + ((constant.c ** 2) * self.r)
        num2 = constant.G * self.M
        deno = (constant.c ** 4) * (self.r ** 3)
        return num1 * num2 / deno

    @property
    def christ_sym1_11(self):
        num = constant.G * self.M
        deno1 = 2 * constant.G * self.M * self.r
        deno2 = (constant.c ** 2) * (self.r ** 2)
        deno = deno1 - deno2
        return num / deno

    @property
    def christ_sym1_22(self):
        return schwarzschild_radius(self.M) - self.r

    @property
    def christ_sym1_33(self):
        return (
            (-1 + (schwarzschild_radius(self.M) / self.r))
            * self.r
            * np.sin(self.theta) ** 2
        )

    @property
    def christ_sym0_01(self):
        num = constant.G * self.M
        deno = self.r * (-2 * constant.G * self.M + (constant.c ** 2) * self.r)
        return num / deno

    @property
    def christ_sym2_21(self):
        return 1 / self.r

    @property
    def christ_sym2_33(self):
        return -1 * np.cos(self.theta) * np.sin(self.theta)

    @property
    def christ_sym3_31(self):
        return 1 / self.r

    @property
    def christ_sym3_32(self):
        return 1 / np.tan(self.theta)

    def get_metric(self):
        """
        Utility to get the Schwarzschild Metric in the form of numpy array.

        Returns
        -------
        metric : ~numpy.array
             Schwarzschild Metric

        """
        metric = np.zeros((4, 4), dtype=float)
        metric[0][0] = self.g00
        metric[1][1] = self.g11
        metric[2][2] = self.g22
        metric[3][3] = self.g33
        return metric
