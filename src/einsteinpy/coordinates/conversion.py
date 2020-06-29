import numpy as np

from einsteinpy.coordinates.utils import (
    bl_to_cartesian_fast,
    cartesian_to_bl_fast,
    cartesian_to_spherical_fast,
    four_position,
    four_velocity,
    spherical_to_cartesian_fast,
    v0,
)
from einsteinpy.metric import BaseMetric


class CartesianConversion:
    """
    Class for conversion to and from Cartesian Coordinates in SI units

    """

    def __init__(self, t, x, y, z, v_x=None, v_y=None, v_z=None):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        x : float
            x-Component of 3-Position
        y : float
            y-Component of 3-Position
        z : float
            z-Component of 3-Position
        v_x : float, optional
            x-Component of 3-Velocity
        v_y : float, optional
            y-Component of 3-Velocity
        v_z : float, optional
            z-Component of 3-Velocity

        """
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        self._velocities_provided = not (
            (v_x is None) or (v_y is None) or (v_z is None)
        )

    def values(self):
        """
        Returns components of the coordinates

        Returns
        -------
        tuple
            4-Tuple containing ``t, x, y, z``

        """
        if self._velocities_provided:
            return self.t, self.x, self.y, self.z, self.v_x, self.v_y, self.v_z
        return self.t, self.x, self.y, self.z

    def convert_spherical(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Spherical Polar Coordinates

        """
        converted = cartesian_to_spherical_fast(
            self.x,
            self.y,
            self.z,
            self.v_x,
            self.v_y,
            self.v_z,
            self._velocities_provided,
        )

        return np.hstack((self.t, converted))

    def convert_bl(self, **kwargs):
        """
        Converts to Boyer-Lindquist Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Boyer-Lindquist Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        alpha = BaseMetric.alpha(M=M, a=a)
        converted = cartesian_to_bl_fast(
            self.x,
            self.y,
            self.z,
            alpha,
            self.v_x,
            self.v_y,
            self.v_z,
            self._velocities_provided,
        )

        return np.hstack((self.t, converted))


class SphericalConversion:
    """
    Class for conversion to and from Spherical Polar Coordinates in SI units

    """

    def __init__(self, t, r, theta, phi, v_r=None, v_th=None, v_p=None):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        r : float
            r-Component of 3-Position
        theta : float
            theta-Component of 3-Position
        phi : float
            phi-Component of 3-Position
        v_r : float, optional
            r-Component of 3-Velocity
        v_th : float, optional
            theta-Component of 3-Velocity
        v_p : float, optional
            phi-Component of 3-Velocity

        """
        self.t = t
        self.r = r
        self.th = theta
        self.p = phi
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        self._velocities_provided = not (
            (v_r is None) or (v_th is None) or (v_p is None)
        )

    def values(self):
        """
        Returns components of the coordinates

        Returns
        -------
        tuple
            4-Tuple containing ``t, r, theta, phi``

        """
        if self._velocities_provided:
            return self.t, self.r, self.th, self.p, self.v_r, self.v_th, self.v_p
        return self.t, self.r, self.th, self.p

    def convert_cartesian(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Other Parameters
        ----------------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Cartesian Coordinates

        """
        converted = spherical_to_cartesian_fast(
            self.r,
            self.th,
            self.p,
            self.v_r,
            self.v_th,
            self.v_p,
            self._velocities_provided,
        )

        return np.hstack((self.t, converted))

    def convert_bl(self, **kwargs):
        """
        Converts to Boyer-Lindquist Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Boyer-Lindquist Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        transformed_cartesian = self.convert_cartesian()
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_bl(M=M, a=a)


class BoyerLindquistConversion:
    """
    Class for conversion to and from Boyer-Lindquist Coordinates in SI units

    """

    def __init__(self, t, r, theta, phi, v_r=None, v_th=None, v_p=None):
        """
        Constructor

        Parameters
        ----------
        t : float
            Time
        r : float
            r-Component of 3-Position
        theta : float
            theta-Component of 3-Position
        phi : float
            phi-Component of 3-Position
        v_r : float, optional
            r-Component of 3-Velocity
        v_th : float, optional
            theta-Component of 3-Velocity
        v_p : float, optional
            phi-Component of 3-Velocity

        """
        self.t = t
        self.r = r
        self.th = theta
        self.p = phi
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        self._velocities_provided = not (
            (v_r is None) or (v_th is None) or (v_p is None)
        )

    def values(self):
        """
        Returns components of the coordinates

        Returns
        -------
        tuple
            4-Tuple containing ``t, r, theta, phi``

        """
        if self._velocities_provided:
            return self.t, self.r, self.th, self.p, self.v_r, self.v_th, self.v_p
        return self.t, self.r, self.th, self.p

    def convert_cartesian(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Cartesian Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        alpha = BaseMetric.alpha(M=M, a=a)
        converted = bl_to_cartesian_fast(
            self.r,
            self.th,
            self.p,
            alpha,
            self.v_r,
            self.v_th,
            self.v_p,
            self._velocities_provided,
        )

        return np.hstack((self.t, converted))

    def convert_spherical(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments
            Expects two arguments, ``M and ``a``, as described below

        Other Parameters
        ----------------
        M : float
            Mass of the gravitating body, \
            around which, spacetime has been defined
        a : float
            Spin Parameter of the gravitating body, \
            around which, spacetime has been defined

        Returns
        -------
        ~numpy.ndarray
            A Length-7 Numpy array, containing the \
            components in Spherical Polar Coordinates

        """
        M, a = kwargs["M"], kwargs["a"]
        transformed_cartesian = self.convert_cartesian(M=M, a=a)
        cart = CartesianConversion(*transformed_cartesian)
        return cart.convert_spherical()
