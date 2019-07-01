"""Important Bodies.
Contains some predefined bodies of the Solar System:
* Sun (☉)
* Earth (♁)
* Moon (☾)
* Mercury (☿)
* Venus (♀)
* Mars (♂)
* Jupiter (♃)
* Saturn (♄)
* Uranus (⛢)
* Neptune (♆)
* Pluto (♇)
and a way to define new bodies (:py:class:`~Body` class).
Data references can be found in :py:mod:`~einsteinpy.constant`
"""

import astropy.units as u

from einsteinpy import constant
from einsteinpy.coordinates import CartesianDifferential


class Body:
    """
    Class to create a generic Body
    """

    @u.quantity_input(mass=u.kg, R=u.km)
    def __init__(
        self,
        name="Generic Body",
        mass=0 * u.kg,
        R=0 * u.km,
        differential=None,
        a=0 * u.m,
        q=0 * u.C,
        parent=None,
    ):
        """
        Parameters
        ----------
        name : str
            Name/ID of the body
        mass : ~astropy.units.kg
            Mass of the body
        R : ~astropy.units
            Radius of the body
        differential : ~einsteinpy.coordinates, optional
            Complete coordinates of the body
        a : float, optional
            Spin factor of massive body
        q : ~astropy.units.C, optional
            Charge on the massive body
        is_attractor : Bool, optional
            To denote is this body is acting as attractor or not
        parent : Body, optional
            The parent object of the body.
        """
        if differential:
            if differential.system == "Cartesian":
                self.pos_vec = [differential.x, differential.y, differential.z]
                self.vel_vec = [differential.v_x, differential.v_y, differential.v_z]
            else:
                self.pos_vec = [differential.r, differential.theta, differential.phi]
                self.vel_vec = [differential.v_r, differential.v_t, differential.v_p]
        self.a = a
        self.R = R
        self.q = q
        self.mass = mass
        self.name = name
        self.coordinates = differential
        self.parent = parent


class _Sun(Body):
    def __init__(self):
        parent = None
        name = "Sun"
        R = constant.R_sun
        mass = constant.Solar_Mass
        super(_Sun, self).__init__(name=name, mass=mass, R=R, parent=parent)


Sun = _Sun()


class _Earth(Body):
    def __init__(self):
        parent = Sun
        name = "Earth"
        R = 6731 * u.km
        mass = 5.97219e24 * u.kg
        super(_Earth, self).__init__(name=name, mass=mass, R=R, parent=parent)


Earth = _Earth()
