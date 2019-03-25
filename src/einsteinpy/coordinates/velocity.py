import numpy as np

import astropy.units as u

class CartesianVelocity:
    """
    Class for representing velocity in Cartesian Coordinates and related transformations.
    """

    @u.quantity_input(x=u.m/u.s, y=u.m/u.s, z=u.m/u.s)
    def __init__(self, vx, vy, vz):
        """
        Constructor.
        Parameters
        ----------
        vx : ~astropy.units
        vy : ~astropy.units
        vz : ~astropy.units
        """
        self.vx = vx
        self.vy = vy
        self.vz = vz

    def to_spherical(self, pos):
        """
        Method for conversion to velocity in spherical coordinates.
        
        Returns
        -------
        spherical : ~einsteinpy.coordinates.core.SphericalVelocity
            Spherical representation of the velocity in Cartesian Coordinates.
        """
        tempvar1 = pos.x ** 2 + pos.y ** 2
        tempvar2 = tempvar1 + pos.z ** 2
        
        vr = (self.vx * pos.x + self.vy * pos.y + self.vz * pos.z) / np.sqrt(tempvar2)
        vtheta = pos.z * (pos.x * self.vx + pos.y * self.vy) - tempvar1 * self.vz
        vphi = -1 * (self.vx * pos.y - self.vy * pos.x) / tempvar1

        return SphericalVelocity(vr, vtheta, vphi)

    def to_bl(self, pos, a):
        """
        Method for conversion of velocity into boyer-lindquist coordinates.
        Returns
        -------
        bl : ~einsteinpy.coordinates.core.BoyerLindquist
            BL representation of the Cartesian Coordinates.
         """
        pos_vec = [pos.x, pos.y, pos.z]
        vel_vec = [self.vx, self.vy, self.vz]
        
        w = float(np.sum(np.square(pos_vec[:3])) - (a ** 2))
        dw_dt = float(2 * np.sum(np.multiply(pos_vec, vel_vec)))
        r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2)))))
        
        vr = (1 / (2 * r)) * (
        (dw_dt / 2)
        +   (
            (w * dw_dt + 4 * (a ** 2) * pos_vec[2] * vel_vec[2])
            / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2))))
            )
        )
        vtheta = (-1 / np.sqrt(1 - np.square(pos_vec[2] / r))) * (
            (vel_vec[2] * r - v_vec[0] * pos_vec[2]) / (r ** 2)
        )
        vphi = (1 / (1 + np.square(pos_vec[1] / pos_vec[0]))) * (
            (vel_vec[1] * pos_vec[0] - vel_vec[0] * pos_vec[1]) / (pos_vec[0] ** 2)
        )

        return BoyerLindquistVelocity(vr, vtheta, vphi)


class SphericalVelocity:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.m/u.s, theta=u.rad/u.s, phi=u.rad/u.sec)
    def __init__(self, vr, vtheta, vphi):
        self.vr = vr
        self.vtheta = vtheta
        self.vphi = vphi

class BoyerLindquistVelocity:
    """
    Class for Spherical Coordinates and related transformations.
    """

    @u.quantity_input(r=u.m/u.s, theta=u.rad/u.s, phi=u.rad/u.s)
    def __init__(self, vr, vtheta, vphi):
        self.vr = vr
        self.vtheta = vtheta
        self.vphi = vphi
