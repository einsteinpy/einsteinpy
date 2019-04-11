import astropy.units as u
import numpy as np

from einsteinpy.coordinates import BoyerLindquist, Cartesian, Spherical


class CartesianDifferential(Cartesian):
    """
    Class for calculating and transforming the velocity
    """

    @u.quantity_input(
        x=u.km, y=u.km, z=u.km, v_x=u.km / u.s, v_y=u.km / u.s, v_z=u.km / u.s
    )
    def __init__(self, x, y, z, v_x, v_y, v_z):
        """
        Constructor.

        Parameters
        ----------
        x : ~astropy.units
        y : ~astropy.units
        z : ~astropy.units
        v_x : ~astropy.units
        v_y : ~astropy.units
        v_z : ~astropy.units
        """
        super(CartesianDifferential, self).__init__(x, y, z)
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z

    def spherical_differential(self):
        """
        Function to convert velocity to spherical coordinates

        Returns
        -------
        SphericalDifferential : ~einsteinpy.coordinates.velocity.SphericalDifferential
            Spherical representation of the velocity in Cartesian Coordinates.

        """
        transformed_spherical = self.to_spherical()
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
            transformed_spherical.r,
            transformed_spherical.theta,
            transformed_spherical.phi,
            v_r,
            v_t,
            v_p,
        )

    def bl_differential(self, a):
        """
        Function to convert velocity to Boyer-Lindquist coordinates

        Returns
        -------
        BoyerLindquistDifferential : ~einsteinpy.coordinates.velocity.BoyerLindquistDifferential
            Boyer-Lindquist representation of the velocity in Cartesian Coordinates.

        """
        transformed_bl = self.to_bl(a)
        w = (self.x ** 2 + self.y ** 2 + self.z ** 2) - (a ** 2)
        dw_dt = 2 * (self.x * self.v_x + self.y * self.v_y + self.z * self.v_z)
        v_r = (1 / (2 * transformed_bl.r)) * (
            (dw_dt / 2)
            + (
                (w * dw_dt + 4 * (a ** 2) * self.z * self.v_z)
                / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2))))
            )
        )
        v_t = (
            (-1 / np.sqrt(1 - np.square(self.z / transformed_bl.r)))
            * ((self.v_z * transformed_bl.r - v_r * self.z) / (transformed_bl.r ** 2))
            * u.rad
        )
        v_p = (
            (1 / (1 + np.square(self.y / self.x)))
            * ((self.v_y * self.x - self.v_x * self.y) / (self.x ** 2))
            * u.rad
        )
        return BoyerLindquistDifferential(
            transformed_bl.r, transformed_bl.theta, transformed_bl.phi, v_r, v_t, v_p
        )


class SphericalDifferential(Spherical):
    """
    Class for calculating and transforming the velocity
    """

    @u.quantity_input(
        r=u.km, theta=u.rad, phi=u.rad, v_r=u.km / u.s, v_t=u.rad / u.s, v_p=u.rad / u.s
    )
    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units
        theta : ~astropy.units
        phi : ~astropy.units
        v_r : ~astropy.units
        v_t : ~astropy.units
        v_p : ~astropy.units
        """
        super(SphericalDifferential, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p

    def cartesian_differential(self):
        """
        Function to convert velocity to cartesian coordinates

        Returns
        -------
        CartesianDifferential : ~einsteinpy.coordinates.velocity.CartesianDifferential
            Cartesian representation of the velocity in Spherical Coordinates.

        """
        transformed_cartesian = self.to_cartesian()
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
        return CartesianDifferential(
            transformed_cartesian.x,
            transformed_cartesian.y,
            transformed_cartesian.z,
            v_x,
            v_y,
            v_z,
        )

    def bl_differential(self, a):
        """
        Function to convert velocity to Boyer-Lindquist coordinates

        Returns
        -------
        BoyerLindquistDifferential : ~einsteinpy.coordinates.velocity.BoyerLindquistDifferential
            Boyer-Lindquist representation of the velocity in Spherical Coordinates.

        """
        transformed_cartesian = self.cartesian_differential()
        return transformed_cartesian.bl_differential(a)


class BoyerLindquistDifferential(BoyerLindquist):
    """
    Class for calculating and transforming the velocity
    """

    @u.quantity_input(
        r=u.km, theta=u.rad, phi=u.rad, v_r=u.km / u.s, v_t=u.rad / u.s, v_p=u.rad / u.s
    )
    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units
        theta : ~astropy.units
        phi : ~astropy.units
        v_r : ~astropy.units
        v_t : ~astropy.units
        v_p : ~astropy.units
        """
        super(BoyerLindquistDifferential, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p

    @u.quantity_input(a=u.km)
    def cartesian_differential(self, a):
        """
        Function to convert velocity to cartesian coordinates

        Parameters
        ----------
        a : float
             a = J/M , the angular momentum per unit mass of the black hole.

        Returns
        -------
        CartesianDifferential : ~einsteinpy.coordinates.velocity.CartesianDifferential
            Cartesian representation of the velocity in Boyer-Lindquist Coordinates.

        """
        transformed_cartesian = self.to_cartesian(a)
        xa = np.sqrt(self.r ** 2 + a ** 2)
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
        return CartesianDifferential(
            transformed_cartesian.x,
            transformed_cartesian.y,
            transformed_cartesian.z,
            v_x,
            v_y,
            v_z,
        )

    @u.quantity_input(a=u.km)
    def spherical_differential(self, a):
        """
        Function to convert velocity to spherical coordinates

        Parameters
        ----------
        a : float
             a = J/M , the angular momentum per unit mass of the black hole.

        Returns
        -------
        SphericalDifferential : ~einsteinpy.coordinates.velocity.SphericalDifferential
            Spherical representation of the velocity in Boyer-Lindquist Coordinates.

        """
        transformed_cartesian = self.cartesian_differential(a)
        return transformed_cartesian.spherical_differential()
