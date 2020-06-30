import astropy.units as u
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

    @u.quantity_input(
        t=u.s, x=u.m, y=u.m, z=u.m, v_x=u.m / u.s, v_y=u.m / u.s, v_z=u.m / u.s
    )
    def __init__(self, t, x, y, z, v_x, v_y, v_z):
        """
        Constructor

        Parameters
        ----------
        t : ~astropy.units.quantity.Quantity
            Time
        x : ~astropy.units.quantity.Quantity
            x-Component of 3-Position
        y : ~astropy.units.quantity.Quantity
            y-Component of 3-Position
        z : ~astropy.units.quantity.Quantity
            z-Component of 3-Position
        v_x : ~astropy.units.quantity.Quantity, optional
            x-Component of 3-Velocity
        v_y : ~astropy.units.quantity.Quantity, optional
            y-Component of 3-Velocity
        v_z : ~astropy.units.quantity.Quantity, optional
            z-Component of 3-Velocity

        """
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        super().__init__(t, x.value, y.value, z.value, v_x.value, v_y.value, v_z.value)
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
        x4 = four_position(self.t.value, self.x.value, self.y.value, self.z.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_x.value, self.v_y.value, self.v_z.value, time_like
        )

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
        x4 = four_position(self.t.value, self.x.value, self.y.value, self.z.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_x.value, self.v_y.value, self.v_z.value, time_like
        )

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
        return SphericalDifferential(
            t * u.s,
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_th * u.rad / u.s,
            v_p * u.rad / u.s,
        )

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
        return BoyerLindquistDifferential(
            t * u.s,
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_th * u.rad / u.s,
            v_p * u.rad / u.s,
        )


class SphericalDifferential(SphericalConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Spherical Polar Coordinates \
    using SI units

    """

    @u.quantity_input(
        t=u.s,
        r=u.m,
        theta=u.rad,
        phi=u.rad,
        v_r=u.m / u.s,
        v_th=u.rad / u.s,
        v_p=u.rad / u.s,
    )
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
            t.value, r.value, theta.value, phi.value, v_r.value, v_th.value, v_p.value
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
        print(self.t, self.r, self.theta, self.phi)
        x4 = four_position(self.t.value, self.r.value, self.theta.value, self.phi.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_r.value, self.v_th.value, self.v_p.value, time_like
        )

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
        print(self.t, self.r, self.theta, self.phi)
        x4 = four_position(self.t.value, self.r.value, self.theta.value, self.phi.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_r.value, self.v_th.value, self.v_p.value, time_like
        )

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
        return CartesianDifferential(
            t * u.s,
            x * u.m,
            y * u.m,
            z * u.m,
            v_x * u.m / u.s,
            v_y * u.m / u.s,
            v_z * u.m / u.s,
        )

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
        return BoyerLindquistDifferential(
            t * u.s,
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_th * u.rad / u.s,
            v_p * u.rad / u.s,
        )


class BoyerLindquistDifferential(BoyerLindquistConversion):
    """
    Class for defining 3-Velocity & 4-Velocity in Boyer-Lindquist Coordinates \
    using SI units

    """

    @u.quantity_input(
        t=u.s,
        r=u.m,
        theta=u.rad,
        phi=u.rad,
        v_r=u.m / u.s,
        v_th=u.rad / u.s,
        v_p=u.rad / u.s,
    )
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
        super().__init__(
            t.value, r.value, theta.value, phi.value, v_r.value, v_th.value, v_p.value
        )
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
        x4 = four_position(self.t.value, self.r.value, self.theta.value, self.phi.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_r.value, self.v_th.value, self.v_p.value, time_like
        )

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
        x4 = four_position(self.t.value, self.r.value, self.theta.value, self.phi.value)
        g_cov_mat = metric.metric_covariant(x4)
        v4 = four_velocity(
            g_cov_mat, self.v_r.value, self.v_th.value, self.v_p.value, time_like
        )

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
        return CartesianDifferential(
            t * u.s,
            x * u.m,
            y * u.m,
            z * u.m,
            v_x * u.m / u.s,
            v_y * u.m / u.s,
            v_z * u.m / u.s,
        )

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
        return SphericalDifferential(
            t * u.s,
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_th * u.rad / u.s,
            v_p * u.rad / u.s,
        )
