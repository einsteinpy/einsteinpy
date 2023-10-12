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

coords: Union[
    CartesianDifferential, SphericalDifferential, BoyerLindquistDifferential, None
] = None
parent: Union["Body", None] = None

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
        return dedent(
            f"""
                 Body(
                    {self.name},
                    {self.parent.name if self.parent else None},
                    {self.mass}, {self.q}, {self.R},
                    {self.coords}
                )"""
        )

    __repr__ = __str__


Sun = Body(name="Sun", mass=Solar_Mass, R=R_sun, parent=None)


Earth = Body(name="Earth", mass=5.97219e24 * u.kg, R=6731 * u.km, parent=Sun)


Moon = Body(name="Moon", mass=7.34767309e22 * u.kg, R=1737.5 * u.km, parent=Earth)


Mercury = Body(name="Mercury", mass=3.285e23 * u.kg, R=2439.7 * u.km, parent=Sun)


Venus = Body(name="Venus", mass=4.867e24 * u.kg, R=6051.8 * u.km, parent=Sun)


Mars = Body(name="Mars", mass=6.39e23 * u.kg, R=3389.5 * u.km, parent=Sun)


Jupiter = Body(name="Jupiter", mass=1.89813e27 * u.kg, R=69911 * u.km, parent=Sun)


Saturn = Body(name="Saturn", mass=5.683e26 * u.kg, R=58232 * u.km, parent=Sun)


Uranus = Body(name="Uranus", mass=8.681e25 * u.kg, R=25362 * u.km, parent=Sun)


Neptune = Body(name="Neptune", mass=1.024e26 * u.kg, R=24622 * u.km, parent=Sun)


Pluto = Body(name="Pluto", mass=1.309e22 * u.kg, R=1183.3 * u.km, parent=Sun)
