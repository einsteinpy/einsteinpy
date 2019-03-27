import astropy.units as u
import numpy as np


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

        sum = x ** 2 + y ** 2 + z ** 2

        distance = sum ** 0.5

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
        w = self.norm() ** 2 - a ** 2
        r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (self.z ** 2)))))
        theta = np.arccos(self.z / r)
        phi = np.arctan2(self.y, self.x)

        return BoyerLindquist(r, theta, phi)


class Spherical:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad)
    def __init__(self, r, theta, phi):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units
        theta : ~astropy.units
        phi : ~astropy.units
        """
        self.r = r
        self.theta = theta
        self.phi = phi

    def to_cartesian(self):
        """
        Method for conversion to cartesian coordinates.

        Returns
        -------
        cartesian : ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the Spherical Coordinates.

        """
        x = self.r * np.cos(self.phi) * np.sin(self.theta)
        y = self.r * np.sin(self.phi) * np.sin(self.theta)
        z = self.r * np.cos(self.theta)

        return Cartesian(x, y, z)

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
            BL representation of the Spherical Coordinates.

        """
        cart = self.to_cartesian()
        bl = cart.to_bl(a)
        return bl


class BoyerLindquist:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.km, theta=u.rad, phi=u.rad)
    def __init__(self, r, theta, phi):
        """
        Constructor.

        Parameters
        ----------
        r : ~astropy.units
        theta : ~astropy.units
        phi : ~astropy.units
        """
        self.r = r
        self.theta = theta
        self.phi = phi

    @u.quantity_input(a=u.km)
    def to_cartesian(self, a):
        """
        Method for conversion to cartesian coordinates.

        Parameters
        ----------
        a : float
             a = J/M , the angular momentum per unit mass of the black hole.

        Returns
        -------
        bl : ~einsteinpy.coordinates.core.Cartesian
            Cartesian representation of the BL Coordinates.

        """
        sin_norm = np.sqrt(self.r ** 2 + a ** 2) * np.sin(self.theta)
        x = sin_norm * np.cos(self.phi)
        y = sin_norm * np.sin(self.phi)
        z = self.r * np.cos(self.theta)
        return Cartesian(x, y, z)

    @u.quantity_input(a=u.km)
    def to_spherical(self, a):
        """
        Method for conversion to spherical coordinates.

        Parameters
        ----------
        a : float
             a = J/M , the angular momentum per unit mass of the black hole.

        Returns
        -------
        spherical : ~einsteinpy.coordinates.core.Spherical
            Spherical representation of the BL Coordinates.

        """
        cart = self.to_cartesian(a)
        spherical = cart.to_spherical()
        return spherical
