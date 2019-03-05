import warnings

import astropy.units as u
import numpy as np

from einsteinpy import constant, utils
from einsteinpy.integrators import RK45
from einsteinpy.utils import schwarzschild_radius as scr

# paper : Fuerst & Woo, 2004 : https://www.aanda.org/articles/aa/pdf/2004/36/aa0814.pdf
# paper : https://arxiv.org/pdf/1601.02063.pdf
# paper : http://www.math.mcgill.ca/gantumur/math599w17/project-kerr.pdf

_c = constant.c.value
_G = constant.G.value


class KerrNull:
    """
    Class for defining methods for Kerr Null Geodesics (Schwarzschild Null Geodesics when a=0)
    """

    def __init__(self, pos_vec, E, initial_p, freq, time, M, a):
        self.M = M
        self.a = a
        self.E = E
        self.initial_p = initial_p
        self.schwarzschild_r = (
            scr(self.M) * (constant.c ** 2) / (constant.G * self.M)
        )  # in normalized units c,G,M = 1, value would be 2
        self.Lz = self.initial_p[2]
        self.initial_sixvec = np.hstack((time.value, pos_vec, initial_p[:2]))

    @classmethod
    def _classmethod_handler(cls, pos_vec, vel_dir, freq, time, M, a):
        E = 1  # constant.h.value * freq
        total_p = 1  # constant.h.value * freq  # total momentum
        # breaking momentum into components
        cs_vel_dir = utils.BLToCartesian_vel(pos_vec, vel_dir, a)
        cs_pos_vec = utils.BLToCartesian_pos(pos_vec, a)
        initial_p = utils.CartesianToBL_vel(
            cs_pos_vec,
            total_p * (cs_vel_dir / np.sqrt(np.sum(np.square(cs_vel_dir)))),
            a,
        )
        return cls(pos_vec, E, initial_p, freq, time.to(u.s), M.to(u.kg), a)

    @classmethod
    @u.quantity_input(M=u.kg, time=u.s)
    def from_BL(cls, pos_vec, vel_direction, freq, time, M, a):
        """
        Constructor (Initialize using Boyer-Lindquist Coordinates/Spherical Coordinates(when a=0))

        Parameters
        ----------
        pos_vec : ~numpy.array
            3 vector position in BL coordinate system
        vel_direction : ~numpy.array
            3 vector initial direction of EM wave. Vector would be normalized internally.
        freq : float
            Frequency of EM wave. Would determine energy and momentum components internally.
        time : ~astropy.units.s
            Starting time
        M : ~astropy.units.kg
            Mass of black hole
        a : float
            Black Hole Spin Factor. Range : [0,1)
        """
        return cls._classmethod_handler(pos_vec, vel_direction, freq, time, M, a)

    @classmethod
    def from_cartesian(cls, pos_vec, vel_direction, freq, time, M, a):
        # """
        # Constructor (Initialize using Cartesian Coordinates)

        # Parameters
        # ----------
        # pos_vec : ~numpy.array
        #     3 vector position in Cartesian coordinate system
        # vel_direction : ~numpy.array
        #     3 vector initial direction of EM wave. Vector would be normalized internally.
        # freq : float
        #     Frequency of EM wave. Would determine energy and momentum components internally.
        # time : ~astropy.units.s
        #     Starting time
        # M : ~astropy.units.kg
        #     Mass of black hole
        # a : float
        #     Black Hole Spin Factor. Range : [0,1)
        # """
        pass

    def sigma(self, r, theta, a):
        """
        Calculate r^2+a^2*cos^2(theta)

        Parameters
        ----------
        r : float
            Component 'r' in BL coordinates
        theta : float
            Component 'theta' in BL coordinates
        a : float
            Constant

        Returns
        -------
        float
            returns r^2+a^2*cos^2(theta)

        """
        return (r ** 2) + ((a * np.cos(theta)) ** 2)

    def delta(self, r, Rs, a):
        """
        Calculate r^2-Rr+a^2

        Parameters
        ----------
        r : float
            Component 'r' in BL coordinates
        Rs : float
            Schwarzschild Radius
        a : float
            Constant

        Returns
        -------
        float
            returns r^2-Rs*r+a^2

        """
        return (r ** 2) - (Rs * r) + (a ** 2)

    def carters_const(self, theta, p_theta, E, Lz):
        return (p_theta ** 2) + (
            (np.cos(theta) ** 2)
            * (((Lz / np.sin(theta)) ** 2) - (self.a ** 2 * E ** 2))
        )

    def kappa(self, theta, p_theta, E, Lz):
        return (
            self.carters_const(theta, p_theta, E, Lz)
            + (Lz ** 2)
            + (self.a ** 2 * E ** 2)
        )

    def f_vec(self, ld, sixvec):
        sigma = self.sigma(sixvec[1], sixvec[2], self.a)
        delta = self.delta(sixvec[1], self.schwarzschild_r, self.a)
        sd = sigma * delta
        # carters_const = self.carters_const(sixvec[2], sixvec[5], self.E, self.Lz)
        kappa = self.kappa(sixvec[2], sixvec[5], self.E, self.Lz)

        def f_t():
            num = (
                2 * sixvec[1] * ((sixvec[1] ** 2) + (self.a ** 2)) * self.E
                - 2 * self.a * self.Lz
            )
            return self.E + (num / sd)

        def f_r():
            return delta * sixvec[4] / sigma

        def f_theta():
            return sixvec[5] / sigma

        def f_phi():
            num1 = 2 * self.a * sixvec[1] * self.E
            num2 = (sigma - 2 * sixvec[1]) * self.Lz / (np.sin(sixvec[2]) ** 2)
            return (num1 + num2) / sd

        def f_p_r():
            term1 = kappa * (1 - sixvec[1])
            term2 = 2 * sixvec[1] * ((sixvec[1] ** 2) + (self.a ** 2)) * self.E
            term3 = -2 * self.a * self.E * self.Lz
            term4 = 2 * (sixvec[4] ** 2) * (sixvec[1] - 1) / sigma
            return ((term1 + term2 + term3) / sd) - term4

        def f_p_theta():
            term1 = np.sin(sixvec[2]) * np.cos(sixvec[2]) / sigma
            term2 = (self.Lz ** 2) / (np.sin(sixvec[2]) ** 4)
            term3 = (self.a ** 2) * (self.E ** 2)
            return term1 * (term2 - term3)

        return np.array(
            [f_t(), f_r(), f_theta(), f_phi(), f_p_r(), f_p_theta()], dtype=float
        )

    def calculate_trajectory(self, end_lambda, OdeMethodKwargs):
        vec_list, lambda_list = list(), list()
        ODE = RK45(self.f_vec, 0.0, self.initial_sixvec, 1e300, **OdeMethodKwargs)
        while ODE.t < end_lambda:
            vec_list.append(ODE.y)
            lambda_list.append(ODE.t)
            ODE.step()

        return (np.array(lambda_list), np.array(vec_list))
