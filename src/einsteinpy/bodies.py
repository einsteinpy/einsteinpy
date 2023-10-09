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

from dataclasses import dataclass

from astropy import units as u

from einsteinpy.constant import R_sun, Solar_Mass

__all__ = ["Body"]


@dataclass
class Body:
    """
    Class to create a generic Body
    """

    name: str = "Generic Body"
    mass: u.kg = 0 * u.kg
    q: u.C = 0 * u.C
    R: u.km = 0 * u.km
    differential = None
    parent: "Body" = None
    """
    Parameters
    ----------
    name : string
        Name or ID of the body
    mass : ~astropy.units.kg
        Mass of the body
    q : ~astropy.units.C, optional
        Charge on the body
    R : ~astropy.units
        Radius of the body
    differential : ~einsteinpy.coordinates.differential.*, optional
        Complete coordinates of the body
    parent : Body, optional
        The parent object of the body
        Useful in case of multibody systems
    """

    def __post_init__(self):
        if self.differential:
            if self.differential.system == "Cartesian":
                self.pos_vec = [
                    self.differential.x,
                    self.differential.y,
                    self.differential.z,
                ]
                self.vel_vec = [
                    self.differential.v_x,
                    self.differential.v_y,
                    self.differential.v_z,
                ]
            else:
                self.pos_vec = [
                    self.differential.r,
                    self.differential.theta,
                    self.differential.phi,
                ]
                self.vel_vec = [
                    self.differential.v_r,
                    self.differential.v_th,
                    self.differential.v_p,
                ]

    def __str__(self):
        return f"Body: ( Name: ({self.name}), Mass: ({self.mass}), Charge: ({self.q}), Radius: ({self.R}), \n \
            Initial Coordinates: ({self.differential}), Parent Body: ({self.parent}) )"

    __repr__ = __str__


@dataclass
class _Sun(Body):
    parent: Body = None
    name: str = "Sun"
    R: u.km = R_sun
    mass: u.kg = Solar_Mass


Sun = _Sun()


@dataclass
class _Earth(Body):
    parent: Body = Sun
    name: str = "Earth"
    R: u.km = 6731 * u.km
    mass: u.kg = 5.97219e24 * u.kg


Earth = _Earth()


class _Moon(Body):
    parent: Body = Earth
    name: str = "Moon"
    R: u.km = 1737.5 * u.km
    mass: u.kg = 7.34767309e22 * u.kg


Moon = _Moon()


class _Mercury(Body):
    parent: Body = Sun
    name: str = "Mercury"
    R: u.km = 2439.7 * u.km
    mass: u.kg = 3.285e23 * u.kg


Mercury = _Mercury()


class _Venus(Body):
    parent: Body = Sun
    name: str = "Venus"
    R: u.km = 6051.8 * u.km
    mass: u.kg = 4.867e24 * u.kg


Venus = _Venus()


class _Mars(Body):
    parent: Body = Sun
    name: str = "Mars"
    R: u.km = 3389.5 * u.km
    mass: u.kg = 6.39e23 * u.kg


Mars = _Mars()


class _Jupiter(Body):
    parent: Body = Sun
    name: str = "Jupiter"
    R: u.km = 69911 * u.km
    mass: u.kg = 1.89813e27 * u.kg


Jupiter = _Jupiter()


class _Saturn(Body):
    parent: Body = Sun
    name: str = "Saturn"
    R: u.km = 58232 * u.km
    mass: u.kg = 5.683e26 * u.kg


Saturn = _Saturn()


class _Uranus(Body):
    parent: Body = Sun
    name: str = "Uranus"
    R: u.km = 25362 * u.km
    mass: u.kg = 8.681e25 * u.kg


Uranus = _Uranus()


class _Neptune(Body):
    parent: Body = Sun
    name: str = "Neptune"
    R: u.km = 24622 * u.km
    mass: u.kg = 1.024e26 * u.kg


Neptune = _Neptune()


class _Pluto(Body):
    parent: Body = Sun
    name: str = "Pluto"
    R: u.km = 1183.3 * u.km
    mass: u.kg = 1.309e22 * u.kg


Pluto = _Pluto()
