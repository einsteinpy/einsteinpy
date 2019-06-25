import astropy.units as u
import numpy as np

from .core import BoyerLindquist, Cartesian, Spherical


class CartesianDifferential(Cartesian):
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
        super(CartesianDifferential, self).__init__(x, y, z)
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z

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
            self.v_x.to(u.m / u.s),
            self.v_y.to(u.m / u.s),
            self.v_z.to(u.m / u.s),
        ]
        return np.hstack(
            (
                super(CartesianDifferential, self).si_values(),
                [e.value for e in element_list],
            )
        )

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
        sph_pos = self.to_spherical()
        n1 = self.x ** 2 + self.y ** 2
        n2 = n1 + self.z ** 2
        v_r = (self.x * self.v_x + self.y * self.v_y + self.z * self.v_z) / np.sqrt(n2)

        v_t = (
            (self.z * (self.x * self.v_x + self.y * self.v_y) - n1 * self.v_z)
            / (n2 * np.sqrt(n1))
            * u.rad
        )

        v_p = (-1 * (self.v_x * self.y - self.x * self.v_y) / n1) * u.rad
        return SphericalDifferential(
            sph_pos.r, sph_pos.theta, sph_pos.phi, v_r, v_t, v_p
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
        bl_pos = self.to_bl(a)
        w = (self.x ** 2 + self.y ** 2 + self.z ** 2) - (a ** 2)
        dw_dt = 2 * (self.x * self.v_x + self.y * self.v_y + self.z * self.v_z)
        v_r = (1 / (2 * bl_pos.r)) * (
            (dw_dt / 2)
            + (
                (w * dw_dt + 4 * (a ** 2) * self.z * self.v_z)
                / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2))))
            )
        )
        v_t = (
            (-1 / np.sqrt(1 - np.square(self.z / bl_pos.r)))
            * ((self.v_z * bl_pos.r - v_r * self.z) / (bl_pos.r ** 2))
            * u.rad
        )
        v_p = (
            (1 / (1 + np.square(self.y / self.x)))
            * ((self.v_y * self.x - self.v_x * self.y) / (self.x ** 2))
            * u.rad
        )
        return BoyerLindquistDifferential(
            bl_pos.r, bl_pos.theta, bl_pos.phi, v_r, v_t, v_p, a
        )


class SphericalDifferential(Spherical):
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
        super(SphericalDifferential, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p

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
            self.v_r.to(u.m / u.s),
            self.v_t.to(u.rad / u.s),
            self.v_p.to(u.rad / u.s),
        ]
        return np.hstack(
            (
                super(SphericalDifferential, self).si_values(),
                [e.value for e in element_list],
            )
        )

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
        cart_pos = self.to_cartesian()
        v_x = (
            np.sin(self.theta) * np.cos(self.phi) * self.v_r
            - self.r * np.sin(self.theta) * np.sin(self.phi) * self.v_p / u.rad
            + self.r * np.cos(self.theta) * np.cos(self.phi) * self.v_t / u.rad
        )
        v_y = (
            np.sin(self.theta) * np.sin(self.phi) * self.v_r
            + self.r * np.cos(self.theta) * np.sin(self.phi) * self.v_t / u.rad
            + self.r * np.sin(self.theta) * np.cos(self.phi) * self.v_p / u.rad
        )
        v_z = (
            np.cos(self.theta) * self.v_r
            - self.r * np.sin(self.theta) * self.v_t / u.rad
        )
        return CartesianDifferential(cart_pos.x, cart_pos.y, cart_pos.z, v_x, v_y, v_z)

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
        transformed_cartesian = self.cartesian_differential()
        return transformed_cartesian.bl_differential(a)


class BoyerLindquistDifferential(BoyerLindquist):
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
        super(BoyerLindquistDifferential, self).__init__(r, theta, phi, a)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p

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
            self.v_r.to(u.m / u.s),
            self.v_t.to(u.rad / u.s),
            self.v_p.to(u.rad / u.s),
        ]
        return np.hstack(
            (
                super(BoyerLindquistDifferential, self).si_values(),
                [e.value for e in element_list],
            )
        )

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
            Cartesian representation of the velocity in Boyer-Lindquist Coordinates.

        """
        cart_pos = self.to_cartesian()
        xa = np.sqrt(self.r ** 2 + self.a ** 2)
        v_x = (
            (self.r * self.v_r * np.sin(self.theta) * np.cos(self.phi) / xa)
            + (xa * np.cos(self.theta) * np.cos(self.phi) * self.v_t / u.rad)
            - (xa * np.sin(self.theta) * np.sin(self.phi) * self.v_p / u.rad)
        )
        v_y = (
            (self.r * self.v_r * np.sin(self.theta) * np.sin(self.phi) / xa)
            + (xa * np.cos(self.theta) * np.sin(self.phi) * self.v_t / u.rad)
            + (xa * np.sin(self.theta) * np.cos(self.phi) * self.v_p / u.rad)
        )
        v_z = (self.v_r * np.cos(self.theta)) - (
            self.r * np.sin(self.theta) * self.v_t / u.rad
        )
        return CartesianDifferential(cart_pos.x, cart_pos.y, cart_pos.z, v_x, v_y, v_z)

    def spherical_differential(self):
        """
        Function to convert velocity to spherical coordinates

        Returns
        -------
        ~einsteinpy.coordinates.velocity.SphericalDifferential
            Spherical representation of the velocity in Boyer-Lindquist Coordinates.

        """
        transformed_cartesian = self.cartesian_differential()
        return transformed_cartesian.spherical_differential()
