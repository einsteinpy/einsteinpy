import astropy.units as u
import numpy as np

from einsteinpy.coordinates.conversion import (
    BoyerLindquistConversion,
    CartesianConversion,
    SphericalConversion,
)


class CartesianDifferential(CartesianConversion):
    """
    Class for calculating and transforming the velocity in Cartesian coordinates.
    """

    @u.quantity_input(
        x=u.km, y=u.km, z=u.km, v_x=u.km / u.s, v_y=u.km / u.s, v_z=u.km / u.s
    )
    def __init__(self, x, y, z, v_x, v_y, v_z):
        """
        Constructor.

        Parameters
        ----------
        x : ~astropy.units.quantity.Quantity
        y : ~astropy.units.quantity.Quantity
        z : ~astropy.units.quantity.Quantity
        v_x : ~astropy.units.quantity.Quantity
        v_y : ~astropy.units.quantity.Quantity
        v_z : ~astropy.units.quantity.Quantity

        """
        self.x = x
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        super().__init__(
            x.si.value, y.si.value, z.si.value, v_x.si.value, v_y.si.value, v_z.si.value
        )
        self.system = "Cartesian"

    def __repr__(self):
        return "Cartesian x: {}, y: {}, z: {}\n" "vx: {}, vy: {}, vz: {}".format(
            self.x, self.y, self.z, self.v_x, self.v_y, self.v_z
        )

    def __str__(self):
        return self.__repr__()

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, m, m, m/s, m/s, m/s)

        """
        element_list = [
            self.x.to(u.m),
            self.y.to(u.m),
            self.z.to(u.m),
            self.v_x.to(u.m / u.s),
            self.v_y.to(u.m / u.s),
            self.v_z.to(u.m / u.s),
        ]
        return np.array([e.value for e in element_list], dtype=float)

    def velocities(self, return_np=False):

        """
        Function for returning velocity.

        Parameters
        ----------
        return_np : bool
            True for numpy array with SI values, False for list with astropy units.
            Defaults to False

        Returns
        -------
        ~numpy.ndarray or list
            Array or list containing velocity.

        """
        if return_np:
            return self.si_values()[3:]
        return [self.v_x, self.v_y, self.v_z]

    def spherical_differential(self):
        """
        Function to convert velocity to spherical coordinates velocity

        Returns
        -------
        ~einsteinpy.coordinates.velocity.SphericalDifferential
            Spherical representation of the velocity in Cartesian Coordinates.

        """
        r, theta, phi, v_r, v_t, v_p = self.convert_spherical()
        return SphericalDifferential(
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_t * u.rad / u.s,
            v_p * u.rad / u.s,
        )

    @u.quantity_input(a=u.km)
    def bl_differential(self, a):
        """
        Function to convert velocity to Boyer-Lindquist coordinates

        Parameters
        ----------
        a : ~astropy.units.quantity.Quantity
            a = J/Mc , the angular momentum per unit mass of the black hole per speed of light.

        Returns
        -------
        ~einsteinpy.coordinates.velocity.BoyerLindquistDifferential
            Boyer-Lindquist representation of the velocity in Cartesian Coordinates.

        """
        r, theta, phi, v_r, v_t, v_p, a = self.convert_bl(a.si.value)
        return BoyerLindquistDifferential(
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_t * u.rad / u.s,
            v_p * u.rad / u.s,
            a * u.m,
        )


class SphericalDifferential(SphericalConversion):
    """
    Class for calculating and transforming the velocity in Spherical coordinates.
    """

    @u.quantity_input(
        r=u.km, theta=u.rad, phi=u.rad, v_r=u.km / u.s, v_t=u.rad / u.s, v_p=u.rad / u.s
    )
    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units.quantity.Quantity
        theta : ~astropy.units.quantity.Quantity
        phi : ~astropy.units.quantity.Quantity
        v_r : ~astropy.units.quantity.Quantity
        v_t : ~astropy.units.quantity.Quantity
        v_p : ~astropy.units.quantity.Quantity

        """
        self.r = r
        self.theta = theta
        self.phi = phi
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p
        super().__init__(
            r.si.value,
            theta.si.value,
            phi.si.value,
            v_r.si.value,
            v_t.si.value,
            v_p.si.value,
        )
        self.system = "Spherical"

    def __repr__(self):
        return "Spherical r: {}, theta: {}, phi: {}\n" "vr: {}, vt: {}, vp: {}".format(
            self.r, self.theta, self.phi, self.v_r, self.v_t, self.v_p
        )

    def __str__(self):
        return self.__repr__()

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, rad, rad, m/s, rad/s, rad/s)

        """
        element_list = [
            self.r.to(u.m),
            self.theta.to(u.rad),
            self.phi.to(u.rad),
            self.v_r.to(u.m / u.s),
            self.v_t.to(u.rad / u.s),
            self.v_p.to(u.rad / u.s),
        ]
        return np.array([e.value for e in element_list], dtype=float)

    def velocities(self, return_np=False):

        """
        Function for returning velocity.

        Parameters
        ----------
        return_np : bool
            True for numpy array with SI values, False for list with astropy units.
            Defaults to False

        Returns
        -------
        ~numpy.ndarray or list
            Array or list containing velocity.

        """

        if return_np:
            return self.si_values()[3:]
        return [self.v_r, self.v_t, self.v_p]

    def cartesian_differential(self):
        """
        Function to convert velocity to cartesian coordinates

        Returns
        -------
        ~einsteinpy.coordinates.velocity.CartesianDifferential
            Cartesian representation of the velocity in Spherical Coordinates.

        """
        x, y, z, v_x, v_y, v_z = self.convert_cartesian()
        return CartesianDifferential(
            x * u.m, y * u.m, z * u.m, v_x * u.m / u.s, v_y * u.m / u.s, v_z * u.m / u.s
        )

    @u.quantity_input(a=u.km)
    def bl_differential(self, a):
        """
        Function to convert velocity to Boyer-Lindquist coordinates

        Parameters
        ----------
        a : ~astropy.units.quantity.Quantity
            a = J/Mc , the angular momentum per unit mass of the black hole per speed of light.

        Returns
        -------
        ~einsteinpy.coordinates.velocity.BoyerLindquistDifferential
            Boyer-Lindquist representation of the velocity in Spherical Coordinates.

        """
        r, theta, phi, v_r, v_t, v_p, a = self.convert_bl(a.si.value)
        return BoyerLindquistDifferential(
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_t * u.rad / u.s,
            v_p * u.rad / u.s,
            a * u.m,
        )


class BoyerLindquistDifferential(BoyerLindquistConversion):
    """
    Class for calculating and transforming the velocity in Boyer-Lindquist coordinates
    """

    @u.quantity_input(
        r=u.km,
        theta=u.rad,
        phi=u.rad,
        v_r=u.km / u.s,
        v_t=u.rad / u.s,
        v_p=u.rad / u.s,
        a=u.km,
    )
    def __init__(self, r, theta, phi, v_r, v_t, v_p, a):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units.quantity.Quantity
        theta : ~astropy.units.quantity.Quantity
        phi : ~astropy.units.quantity.Quantity
        v_r : ~astropy.units.quantity.Quantity
        v_t : ~astropy.units.quantity.Quantity
        v_p : ~astropy.units.quantity.Quantity
        a : ~astropy.units.quantity.Quantity

        """
        self.r = r
        self.theta = theta
        self.phi = phi
        self.a = a
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p
        super().__init__(
            r.si.value,
            theta.si.value,
            phi.si.value,
            v_r.si.value,
            v_t.si.value,
            v_p.si.value,
            a.si.value,
        )
        self.system = "BoyerLindquist"

    def __repr__(self):
        return (
            "Boyer-Lindquist r: {}, theta: {}, phi: {}\n"
            "vr: {}, vt: {}, vp: {}\n"
            "a: {}".format(
                self.r, self.theta, self.phi, self.v_r, self.v_t, self.v_p, self.a
            )
        )

    def __str__(self):
        return self.__repr__()

    def si_values(self):
        """
        Function for returning values in SI units.

        Returns
        -------
        ~numpy.ndarray
            Array containing values in SI units (m, rad, rad, m/s, rad/s, rad/s)

        """
        element_list = [
            self.r.to(u.m),
            self.theta.to(u.rad),
            self.phi.to(u.rad),
            self.v_r.to(u.m / u.s),
            self.v_t.to(u.rad / u.s),
            self.v_p.to(u.rad / u.s),
        ]
        return np.array([e.value for e in element_list], dtype=float)

    def velocities(self, return_np=False):

        """
        Function for returning velocity.

        Parameters
        ----------
        return_np : bool
            True for numpy array with SI values, False for list with astropy units.
            Defaults to False

        Returns
        -------
        ~numpy.ndarray or list
            Array or list containing velocity.

        """
        if return_np:
            return self.si_values()[3:6]
        return [self.v_r, self.v_t, self.v_p]

    def cartesian_differential(self):
        """
        Function to convert velocity to cartesian coordinates

        Returns
        -------
        ~einsteinpy.coordinates.velocity.CartesianDifferential
            Cartesian representation of the velocity in Boyer-Lindquist Coordinates.

        """
        x, y, z, v_x, v_y, v_z = self.convert_cartesian()
        return CartesianDifferential(
            x * u.m, y * u.m, z * u.m, v_x * u.m / u.s, v_y * u.m / u.s, v_z * u.m / u.s
        )

    def spherical_differential(self):
        """
        Function to convert velocity to spherical coordinates

        Returns
        -------
        ~einsteinpy.coordinates.velocity.SphericalDifferential
            Spherical representation of the velocity in Boyer-Lindquist Coordinates.

        """
        r, theta, phi, v_r, v_t, v_p = self.convert_spherical()
        return SphericalDifferential(
            r * u.m,
            theta * u.rad,
            phi * u.rad,
            v_r * u.m / u.s,
            v_t * u.rad / u.s,
            v_p * u.rad / u.s,
        )
