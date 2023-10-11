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
from textwrap import dedent
from typing import Union

from astropy import units as u

from einsteinpy.constant import R_sun, Solar_Mass
from einsteinpy.coordinates.differential import (
    BoyerLindquistDifferential,
    CartesianDifferential,
    SphericalDifferential,
)

__all__ = ["Body"]


@dataclass
class Body:
    """
    Class to create a generic Body

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
    coords : ~einsteinpy.coordinates.differential.*, optional
        Complete coordinates of the body
    parent : Body, optional
        The parent object of the body
        Useful in case of multibody systems
    """

    name: str = "Generic Body"
    mass: u.kg = 0 * u.kg
    q: u.C = 0 * u.C
    R: u.km = 0 * u.km
    coords = None
    parent: "Body" = None

    def __post_init__(self):
        @u.quantity_input(mass=u.kg, q=u.C, R=u.km)
        def check_units(mass, q, R):
            return mass, q, R

        check_units(self.mass, self.q, self.R)
        diff = self.coords
        if diff:
            if diff.system == "Cartesian":
                self.pos_vec = [
                    diff.x,
                    diff.y,
                    diff.z,
                ]
                self.vel_vec = [
                    diff.v_x,
                    diff.v_y,
                    diff.v_z,
                ]
            else:
                self.pos_vec = [
                    diff.r,
                    diff.theta,
                    diff.phi,
                ]
                self.vel_vec = [
                    diff.v_r,
                    diff.v_th,
                    diff.v_p,
                ]

    def __str__(self):
        return dedent(f"""
                 Body(
                    {self.name},
                    {self.parent.name if self.parent else None},
                    {self.mass}, {self.q}, {self.R},
                    {self.coords}
                )"""
            )

    __repr__ = __str__


@dataclass
class _Sun(Body):
    parent: Body = None
    name: str = "Sun"
    R: u.km = R_sun
    mass: u.kg = Solar_Mass


Sun = Body(name="Sun", mass=Solar_Mass, R=R_sun, parent=None)


@dataclass
class _Earth(Body):
    parent: Body = Sun
    name: str = "Earth"
    R: u.km = 6731 * u.km
    mass: u.kg = 5.97219e24 * u.kg


Earth = Body(name="Earth", mass=5.97219e24 * u.kg, R=6731 * u.km, parent=Sun)


class _Moon(Body):
    parent: Body = Earth
    name: str = "Moon"
    R: u.km = 1737.5 * u.km
    mass: u.kg = 7.34767309e22 * u.kg


Moon = Body(name="Moon", mass=7.34767309e22 * u.kg, R=1737.5 * u.km, parent=Earth)


class _Mercury(Body):
    parent: Body = Sun
    name: str = "Mercury"
    R: u.km = 2439.7 * u.km
    mass: u.kg = 3.285e23 * u.kg


Mercury = Body(name="Mercury", mass=3.285e23 * u.kg, R=2439.7 * u.km, parent=Sun)


class _Venus(Body):
    parent: Body = Sun
    name: str = "Venus"
    R: u.km = 6051.8 * u.km
    mass: u.kg = 4.867e24 * u.kg


Venus = Body(name="Venus", mass=4.867e24 * u.kg, R=6051.8 * u.km, parent=Sun)


class _Mars(Body):
    parent: Body = Sun
    name: str = "Mars"
    R: u.km = 3389.5 * u.km
    mass: u.kg = 6.39e23 * u.kg


Mars = Body(name="Mars", mass=6.39e23 * u.kg, R=3389.5 * u.km, parent=Sun)


class _Jupiter(Body):
    parent: Body = Sun
    name: str = "Jupiter"
    R: u.km = 69911 * u.km
    mass: u.kg = 1.89813e27 * u.kg


Jupiter = Body(name="Jupiter", mass=1.89813e27 * u.kg, R=69911 * u.km, parent=Sun)


class _Saturn(Body):
    parent: Body = Sun
    name: str = "Saturn"
    R: u.km = 58232 * u.km
    mass: u.kg = 5.683e26 * u.kg


Saturn = Body(name="Saturn", mass=5.683e26 * u.kg, R=58232 * u.km, parent=Sun)


class _Uranus(Body):
    parent: Body = Sun
    name: str = "Uranus"
    R: u.km = 25362 * u.km
    mass: u.kg = 8.681e25 * u.kg


Uranus = Body(name="Uranus", mass=8.681e25 * u.kg, R=25362 * u.km, parent=Sun)


class _Neptune(Body):
    parent: Body = Sun
    name: str = "Neptune"
    R: u.km = 24622 * u.km
    mass: u.kg = 1.024e26 * u.kg


Neptune = Body(name="Neptune", mass=1.024e26 * u.kg, R=24622 * u.km, parent=Sun)


class _Pluto(Body):
    parent: Body = Sun
    name: str = "Pluto"
    R: u.km = 1183.3 * u.km
    mass: u.kg = 1.309e22 * u.kg


Pluto = Body(name="Pluto", mass=1.309e22 * u.kg, R=1183.3 * u.km, parent=Sun)
