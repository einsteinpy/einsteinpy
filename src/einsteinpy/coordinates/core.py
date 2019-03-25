import numpy as np

import astropy.units as u

class Cartesian:
    """
    Class for Cartesian Coordinates and related transformations.
    """

    @u.quantity_input(x=u.km, y=u.km, z=u.km)
    def __init__(self, x, y, z):
        """
        Constructor.

        Parameters
        ----------
        x : ~astropy.units
        y : ~astropy.units
        z : ~astropy.units
        """
        self.x = x
        self.y = y
        self.z = z

    def norm(self):
        """
        Function for finding euclidean norm of a vector.

        Returns
        -------
        distance : ~astropy.units
            Euclidean norm with units.

        """
        o = Cartesian(0 * u.km, 0 * u.km, 0 * u.km)

        x = self.x - o.x
        y = self.y - o.y
        z = self.z - o.z

        sum = x**2 + y**2 + z**2

        distance = sum**0.5

        return distance

    def dot(self, target):
        """
        Dot product of two vectors.

        Parameters
        ----------
        target: ~einsteipy.coordinates.core.Cartesian

        Returns
        -------
        dotproduct : ~astropy.units
            Dot product with units

        """

        x = self.x * target.x
        y = self.y * target.y
        z = self.z * target.z

        return x + y + z

    def to_spherical(self):
        """
        Method for conversion to spherical coordinates.

        Returns
        -------
        spherical : ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the Cartesian Coordinates.

        """
        r = self.norm()
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return Spherical(r, theta, phi)

    @u.quantity_input(a=u.km)
    def to_bl(self, a):
        """
        Method for conversion to boyer-lindquist coordinates.

        Parameters
        ----------
        a : float
             a = J/M , the angular momentum per unit mass of the black hole.

        Returns
        -------
        bl : ~einsteinpy.coordinates.core.BoyerLindquist
            BL representation of the Cartesian Coordinates.

        """
        w = self.norm()**2 - a ** 2
        r = np.sqrt(
            0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2))))
        )
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return BoyerLindquist(r, theta, phi)


class Spherical:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad)
    def __init__(self, r, theta, phi):
        self.r = r
        self.theta = theta
        self.phi = phi


class BoyerLindquist:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad)
    def __init__(self, r, theta, phi):
        self.r = r
        self.theta = theta
        self.phi = phi

