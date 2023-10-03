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

from astropy import units as u
from einsteinpy import constant

__all__ = ["Body"]

class Body:
    """
    Class to create a generic Body
    """

    @u.quantity_input(mass=u.kg, R=u.km)
    def __init__(
        self,
        name="Generic Body",
        mass=0 * u.kg,
        q=0 * u.C,
        R=0 * u.km,
        differential=None,
        parent=None,
    ):
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
        if differential:
            if differential.system == "Cartesian":
                self.pos_vec = [differential.x, differential.y, differential.z]
                self.vel_vec = [differential.v_x, differential.v_y, differential.v_z]
            else:
                self.pos_vec = [differential.r, differential.theta, differential.phi]
                self.vel_vec = [differential.v_r, differential.v_th, differential.v_p]
        self.name = name
        self.mass = mass
        self.q = q
        self.R = R
        self.coords = differential
        self.parent = parent

    def __str__(self):
        return f"Body: ( Name: ({self.name}), Mass: ({self.mass}), Charge: ({self.q})', Radius: ({self.R}), \n \
            Initial Coordinates: ({self.coords}), Parent Body: ({self.parent}) )"

    def __repr__(self):
        return f"Body: ( Name: ({self.name}), Mass: ({self.mass}), Charge: ({self.q})', Radius: ({self.R}), \n \
            Initial Coordinates: ({self.coords}), Parent Body: ({self.parent}) )"


# Define a dictionary of predefined bodies
predefined_bodies = {
    "Sun": {
        "parent": None,
        "R": constant.R_sun,
        "mass": constant.Solar_Mass,
    },
    "Earth": {
        "parent": "Sun",
        "R": 6731 * u.km,
        "mass": 5.97219e24 * u.kg,
    },
    "Moon": {
        "parent": "Earth",
        "R": 1737.5 * u.km,
        "mass": 7.34767309e22 * u.kg,
    },
    "Mercury": {
        "parent": "Sun",
        "R": 2439.7 * u.km,
        "mass": 3.285e23 * u.kg,
    },
    "Venus": {
        "parent": "Sun",
        "R": 6051.8 * u.km,
        "mass": 4.867e24 * u.kg,
    },
    "Mars": {
        "parent": "Sun",
        "R": 3389.5 * u.km,
        "mass": 6.39e23 * u.kg,
    },
    "Jupiter": {
        "parent": "Sun",
        "R": 69911 * u.km,
        "mass": 1.89813e27 * u.kg,
    },
    "Saturn": {
        "parent": "Sun",
        "R": 58232 * u.km,
        "mass": 5.683e26 * u.kg,
    },
    "Uranus": {
        "parent": "Sun",
        "R": 25362 * u.km,
        "mass": 8.681e25 * u.kg,
    },
    "Neptune": {
        "parent": "Sun",
        "R": 24622 * u.km,
        "mass": 1.024e26 * u.kg,
    },
    "Pluto": {
        "parent": "Sun",
        "R": 1183.3 * u.km,
        "mass": 1.309e22 * u.kg,
    },
}

# Create instances of predefined bodies
for body_name, body_info in predefined_bodies.items():
    body_instance = Body(name=body_name, parent=body_info["parent"], R=body_info["R"], mass=body_info["mass"])
    globals()[body_name] = body_instance
