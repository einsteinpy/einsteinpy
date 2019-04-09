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
        super(CartesianDifferential, self).__init__(x, y, z)
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z

    def spherical_differential(self):
        """
        Function to convert cartesian to spherical coordinates
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

    def boyerlindquistdifferential(self, a):
        """
        Function to convert velocity in Cartesian coordinates to velocity in Boyer-Lindquist coordinates
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
        v_t = (-1 / np.sqrt(1 - np.square(self.z / transformed_bl.r))) * (
            (self.v_z * transformed_bl.r - v_r * self.z) / (transformed_bl.r ** 2)
        )
        v_p = (1 / (1 + np.square(self.y / self.x))) * (
            (self.v_y * self.x - self.v_x * self.y) / (self.x ** 2)
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
        super(SphericalDifferential, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p

    def cartesian_differential(self):
        """
        Function to convert spherical to cartesian coordinates
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

    def boyerlindquistdifferential(self, a):
        """
        Function to convert velocity in spherical coordinates to velocity in Boyer-Lindquist coordinates
        """
        transformed_cartesian = self.cartesian_differential()
        return transformed_cartesian.boyerlindquistdifferential(a)


class BoyerLindquistDifferential(BoyerLindquist):
    """
    Class for calculating and transforming the velocity
    """

    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        super(BoyerLindquistDifferential, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p
