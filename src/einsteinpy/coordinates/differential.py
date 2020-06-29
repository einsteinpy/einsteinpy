import numpy as np

from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)
from einsteinpy.coordinates.utils import four_position, four_velocity


class CartesianDifferential(CartesianConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Cartesian Coordinates \
    using SI units

    """

    def __init__(self, t, x, y, z, v_x, v_y, v_z):
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
        super().__init__(t, x, y, z, v_x, v_y, v_z)
        self.system = "Cartesian"

    def __str__(self):
        return f"Cartesian Coordinates: \n \
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})\n \
            v_x: {self.v_x}, v_y: {self.v_y}, v_z: {self.v_z}"

    def __repr__(self):
        return f"Cartesian Coordinates: \n \
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})\n \
            v_x: {self.v_x}, v_y: {self.v_y}, v_z: {self.v_z}"

    def velocity(self, metric, time_like=True):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray :
            Array, containing Velocity 4-Vector in SI units

        """
        x4 = four_position(self.t, self.x, self.y, self.z)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_x, self.v_y, self.v_z, time_like)

        return v4

    def state(self, metric, time_like=True):
        """
        Returns the State Vector, pertaining to the test particle, \
        whose coordinates have been provided

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        state : ~numpy.ndarray
            The State Vector of the test particle
            Length-8 Array

        """
        x4 = four_position(self.t, self.x, self.y, self.z)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_x, self.v_y, self.v_z, time_like)

        state = np.hstack((x4, v4))

        return state

    def spherical_differential(self, **kwargs):
        """
        Converts to Spherical Polar Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.differential.SphericalDifferential
            Spherical Polar representation of velocity

        """
        t, r, theta, phi, v_r, v_th, v_p = self.convert_spherical()
        return SphericalDifferential(t, r, theta, phi, v_r, v_th, v_p)

    def bl_differential(self, **kwargs):
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
        ~einsteinpy.coordinates.differential.BoyerLindquistDifferential
            Boyer-Lindquist representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        t, r, theta, phi, v_r, v_th, v_p = self.convert_bl(M=M, a=a)
        return BoyerLindquistDifferential(t, r, theta, phi, v_r, v_th, v_p)


class SphericalDifferential(SphericalConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Spherical Polar Coordinates \
    using SI units

    """

    def __init__(self, t, r, theta, phi, v_r, v_th, v_p):
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
        self.theta = theta
        self.phi = phi
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        super().__init__(
            t, r, theta, phi, v_r, v_th, v_p,
        )
        self.system = "Spherical"

    def __str__(self):
        return f"Spherical Polar Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n \
            v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def __repr__(self):
        return f"Spherical Polar Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n \
            v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def velocity(self, metric, time_like=True):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray :
            Array, containing Velocity 4-Vector in SI units

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_r, self.v_th, self.v_p, time_like)

        return v4

    def state(self, metric, time_like=True):
        """
        Returns the State Vector, pertaining to the test particle, \
        whose coordinates have been provided

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        state : ~numpy.ndarray
            The State Vector of the test particle
            Length-8 Array

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_r, self.v_th, self.v_p, time_like)

        state = np.hstack((x4, v4))

        return state

    def cartesian_differential(self, **kwargs):
        """
        Converts to Cartesian Coordinates

        Parameters
        ----------
        **kwargs : dict
            Keyword Arguments

        Returns
        -------
        ~einsteinpy.coordinates.differential.CartesianDifferential
            Cartesian representation of velocity

        """
        t, x, y, z, v_x, v_y, v_z = self.convert_cartesian()
        return CartesianDifferential(t, x, y, z, v_x, v_y, v_z)

    def bl_differential(self, **kwargs):
        """
        Converts to Boyer-Lindquist coordinates

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
        ~einsteinpy.coordinates.differential.BoyerLindquistDifferential
            Boyer-Lindquist representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        t, r, theta, phi, v_r, v_th, v_p = self.convert_bl(M=M, a=a)
        return BoyerLindquistDifferential(t, r, theta, phi, v_r, v_th, v_p)


class BoyerLindquistDifferential(BoyerLindquistConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Boyer-Lindquist Coordinates \
    using SI units

    """

    def __init__(self, t, r, theta, phi, v_r, v_th, v_p):
        """
        Constructor.

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
        self.theta = theta
        self.phi = phi
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        super().__init__(t, r, theta, phi, v_r, v_th, v_p)
        self.system = "BoyerLindquist"

    def __str__(self):
        return f"Boyer-Lindquist Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n \
            v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def __repr__(self):
        return f"Boyer-Lindquist Coordinates: \n \
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n \
            v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def velocity(self, metric, time_like=True):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        ~numpy.ndarray :
            Array, containing Velocity 4-Vector in SI units

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_r, self.v_th, self.v_p, time_like)

        return v4

    def state(self, metric, time_like=True):
        """
        Returns the State Vector, pertaining to the test particle, \
        whose coordinates have been provided

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined
        time_like : bool, optional
            To determine, if the coordinates are for a Time-like or \
            a Null-like Geodesic
            Defaults to ``True``

        Returns
        -------
        state : ~numpy.ndarray
            The State Vector of the test particle
            Length-8 Array

        """
        x4 = four_position(self.t, self.r, self.theta, self.phi)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(g_cov_mat, self.v_r, self.v_th, self.v_p, time_like)

        state = np.hstack((x4, v4))

        return state

    def cartesian_differential(self, **kwargs):
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
        ~einsteinpy.coordinates.differentia.CartesianDifferential
            Cartesian representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        t, x, y, z, v_x, v_y, v_z = self.convert_cartesian(M=M, a=a)
        return CartesianDifferential(t, x, y, z, v_x, v_y, v_z)

    def spherical_differential(self, **kwargs):
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
        ~einsteinpy.coordinates.differentia.SphericalDifferential
            Spherical representation of velocity

        """
        M, a = kwargs["M"], kwargs["a"]
        t, r, theta, phi, v_r, v_th, v_p = self.convert_spherical(M=M, a=a)
        return SphericalDifferential(t, r, theta, phi, v_r, v_th, v_p)
