import numpy as np
from astropy import units as u

from einsteinpy import constant, metric
from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)
from einsteinpy.coordinates.utils import v0
from einsteinpy.utils import CoordinateError

_c = constant.c.value


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
        super().__init__(
            t.si.value,
            x.si.value,
            y.si.value,
            z.si.value,
            v_x.si.value,
            v_y.si.value,
            v_z.si.value,
        )
        self.t = t
        self.x = x
        self.y = y
        self.z = z
        self._v_t = None
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        self.system = "Cartesian"

    def __str__(self):
        return f"Cartesian Coordinates: \n\
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})\n\
            v_t: {self.v_t}, v_x: {self.v_x}, v_y: {self.v_y}, v_z: {self.v_z}"

    def __repr__(self):
        return f"Cartesian Coordinates: \n\
            t = ({self.t}), x = ({self.x}), y = ({self.y}), z = ({self.z})\n\
            v_t: {self.v_t}, v_x: {self.v_x}, v_y: {self.v_y}, v_z: {self.v_z}"

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (_c * self.t.si.value, self.x.si.value, self.y.si.value, self.z.si.value)

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Cartesian Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                f"Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Cartesian Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.v_x.si.value, self.v_y.si.value, self.v_z.si.value)

        self._v_t = v_t * u.m / u.s

    def velocity(self, metric):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined

        Returns
        -------
        tuple
            4-Tuple, containing Velocity 4-Vector in SI units

        """
        # Setting _v_t
        self.v_t = (metric,)

        return (
            self._v_t.value,
            self.v_x.si.value,
            self.v_y.si.value,
            self.v_z.si.value,
        )

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
        super().__init__(
            t.si.value,
            r.si.value,
            theta.si.value,
            phi.si.value,
            v_r.si.value,
            v_th.si.value,
            v_p.si.value,
        )
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
        self._v_t = None
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        self.system = "Spherical"

    def __str__(self):
        return f"Spherical Polar Coordinates: \n\
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n\
            v_t: {self.v_t}, v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def __repr__(self):
        return f"Spherical Polar Coordinates: \n\
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n\
            v_t: {self.v_t}, v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.t.si.value,
            self.r.si.value,
            self.theta.si.value,
            self.phi.si.value,
        )

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Sperical Polar Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                f"Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Spherical Polar Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.v_r.si.value, self.v_th.si.value, self.v_p.si.value)

        self._v_t = v_t * u.m / u.s

    def velocity(self, metric):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined

        Returns
        -------
        tuple
            4-Tuple, containing Velocity 4-Vector in SI units

        """
        # Setting _v_t
        self.v_t = (metric,)

        return (
            self._v_t.value,
            self.v_r.si.value,
            self.v_th.si.value,
            self.v_p.si.value,
        )

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
        super().__init__(
            t.si.value,
            r.si.value,
            theta.si.value,
            phi.si.value,
            v_r.si.value,
            v_th.si.value,
            v_p.si.value,
        )
        self.t = t
        self.r = r
        self.theta = theta
        self.phi = phi
        self._v_t = None
        self.v_r = v_r
        self.v_th = v_th
        self.v_p = v_p
        self.system = "BoyerLindquist"

    def __str__(self):
        return f"Boyer-Lindquist Coordinates: \n\
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n\
            v_t: {self.v_t}, v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def __repr__(self):
        return f"Boyer-Lindquist Coordinates: \n\
            t = ({self.t}), r = ({self.r}), theta = ({self.theta}), phi = ({self.phi})\n\
            v_t: {self.v_t}, v_r: {self.v_r}, v_th: {self.v_th}, v_p: {self.v_p}"

    def position(self):
        """
        Returns Position 4-Vector in SI units

        Returns
        -------
        tuple
            4-Tuple, containing Position 4-Vector in SI units

        """
        return (
            _c * self.t.si.value,
            self.r.si.value,
            self.theta.si.value,
            self.phi.si.value,
        )

    @property
    def v_t(self):
        """
        Returns the Timelike component of 4-Velocity

        """
        return self._v_t

    @v_t.setter
    def v_t(self, args):
        """
        Sets the value of the Time-like component of 4-Velocity

        Parameters
        ----------
        args : tuple
            1-tuple containing the ~einsteinpy.metric.* object, \
            in which the coordinates are defined

        Raises
        ------
        CoordinateError
            If ``metric`` object has been instantiated with a coordinate system, \
            other than Boyer-Lindquist Coordinates.

        """
        g = args[0]
        if self.system != g.coords.system:
            raise CoordinateError(
                "Metric object has been instantiated with a coordinate system, ( {g.coords.system} )"
                " other than Boyer-Lindquist Coordinates."
            )

        g_cov_mat = g.metric_covariant(self.position())

        v_t = v0(g_cov_mat, self.v_r.si.value, self.v_th.si.value, self.v_p.si.value)

        self._v_t = v_t * u.m / u.s

    def velocity(self, metric):
        """
        Returns Velocity 4-Vector in SI units

        Parameters
        ----------
        metric : ~einsteinpy.metric.*
            Metric object, in which the coordinates are defined

        Returns
        -------
        tuple
            4-Tuple, containing Velocity 4-Vector in SI units

        """
        # Setting _v_t
        self.v_t = (metric,)

        return (
            self._v_t.value,
            self.v_r.si.value,
            self.v_th.si.value,
            self.v_p.si.value,
        )

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
