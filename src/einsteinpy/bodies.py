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
        a : ~astropy.units.m, optional
            Spin factor of massive body. Should be less than half of schwarzschild radius.
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

    def __repr__(self):
        return (
            "'Body ( name: ({0}), mass: ({1}), radius: ({2}), coordinates: ({3}), spin factor: ({4}), charge: ({"
            "5}) )'".format(
                self.name, self.mass, self.R, self.coordinates, self.a, self.q
            )
        )

    def __str__(self):
        return (
            "Body ( name: ({0}), mass: ({1}), radius: ({2}), coordinates: ({3}), spin factor: ({4}), charge: ({"
            "5}) )".format(
                self.name, self.mass, self.R, self.coordinates, self.a, self.q
            )
        )


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


class _Moon(Body):
    def __init__(self):
        parent = Earth
        name = "Moon"
        R = 1737.5 * u.km
        mass = 7.34767309e22 * u.kg
        super(_Moon, self).__init__(name=name, mass=mass, R=R, parent=parent)


Moon = _Moon()


class _Mercury(Body):
    def __init__(self):
        parent = Sun
        name = "Mercury"
        R = 2439.7 * u.km
        mass = 3.285e23 * u.kg
        super(_Mercury, self).__init__(name=name, mass=mass, R=R, parent=parent)


Mercury = _Mercury()


class _Venus(Body):
    def __init__(self):
        parent = Sun
        name = "Venus"
        R = 6051.8 * u.km
        mass = 4.867e24 * u.kg
        super(_Venus, self).__init__(name=name, mass=mass, R=R, parent=parent)


Venus = _Venus()


class _Mars(Body):
    def __init__(self):
        parent = Sun
        name = "Mars"
        R = 3389.5 * u.km
        mass = 6.39e23 * u.kg
        super(_Mars, self).__init__(name=name, mass=mass, R=R, parent=parent)


Mars = _Mars()


class _Jupiter(Body):
    def __init__(self):
        parent = Sun
        name = "Jupiter"
        R = 69911 * u.km
        mass = 1.89813e27 * u.kg
        super(_Jupiter, self).__init__(name=name, mass=mass, R=R, parent=parent)


Jupiter = _Jupiter()


class _Saturn(Body):
    def __init__(self):
        parent = Sun
        name = "Saturn"
        R = 58232 * u.km
        mass = 5.683e26 * u.kg
        super(_Saturn, self).__init__(name=name, mass=mass, R=R, parent=parent)


Saturn = _Saturn()


class _Uranus(Body):
    def __init__(self):
        parent = Sun
        name = "Uranus"
        R = 25362 * u.km
        mass = 8.681e25 * u.kg
        super(_Uranus, self).__init__(name=name, mass=mass, R=R, parent=parent)


Uranus = _Uranus()


class _Neptune(Body):
    def __init__(self):
        parent = Sun
        name = "Neptune"
        R = 24622 * u.km
        mass = 1.024e26 * u.kg
        super(_Neptune, self).__init__(name=name, mass=mass, R=R, parent=parent)


Neptune = _Neptune()


class _Pluto(Body):
    def __init__(self):
        parent = Sun
        name = "Pluto"
        R = 1183.3 * u.km
        mass = 1.309e22 * u.kg
        super(_Pluto, self).__init__(name=name, mass=mass, R=R, parent=parent)


Pluto = _Pluto()
