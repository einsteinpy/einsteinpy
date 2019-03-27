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
        super(Cartesian, self).__init__(x, y, z)
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z

    def spherical_differential(self):
        transformed_spherical = self.to_spherical()
        v_r = (
            np.sin(self.y) * np.cos(self.z) * self.v_x
            - self.x * np.sin(self.y) * np.sin(self.z) * self.v_z
            + self.x * np.cos(self.y) * np.cos(self.z) * self.v_y
        )
        v_t = (
            np.sin(self.y) * np.sin(self.z) * self.v_x
            + self.x * np.cos(self.y) * np.sin(self.z) * self.v_y
            + self.x * np.sin(self.y) * np.cos(self.z) * self.v_z
        )
        v_p = np.cos(self.y) * self.v_x - self.x * np.sin(self.y) * self.v_y
        return SphericalDifferential(
            transformed_spherical.x,
            transformed_spherical.y,
            transformed_spherical.z,
            v_r,
            v_t,
            v_p,
        )


class SphericalDifferential(Spherical):
    """
    Class for calculating and transforming the velocity
    """

    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        super(Spherical, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p


class BoyerLindquistDifferential(BoyerLindquist):
    """
    Class for calculating and transforming the velocity
    """

    def __init__(self, r, theta, phi, v_r, v_t, v_p):
        super(BoyerLindquist, self).__init__(r, theta, phi)
        self.v_r = v_r
        self.v_t = v_t
        self.v_p = v_p
