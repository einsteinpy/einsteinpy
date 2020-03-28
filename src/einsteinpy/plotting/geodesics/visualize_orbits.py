from astroquery.jplhorizons import Horizons

from einsteinpy.metric import Schwarzschild
from einsteinpy.utils import schwarzschild_radius


class VisualizeOrbits:
    """
    Class for Visualizing Orbits using JPLHorizons.
    """

    def __init__(self, Id, location, epochs):
        self.Id = Id
        self.location = location
        self.epochs = epochs

    def create_object(self):
        obj = Horizons(id=self.Id, location=self.location, epochs=self.epochs, id_type='id')
        return obj

    def velocity(self):
        # Calculate velocity of the object with respect to sun and the observer.
        # UNITS: Km/s and Km/s
        print(self.create_object().ephemerides(quantities=22))

    def S_distance(self):
        # Calculates "Sub solar" point positional angle and
        # the angular distance from the sub-observer point (center of disk).
        # UNITS: Degrees and Arcseconds
        print(self.create_object().ephemerides(quantities=16))

    def N_distance(self):
        # Calculates North Pole positional angle and the
        # angular distance from the sub-observer point (center of disk).
        # UNITS: Degrees and Arcseconds
        print(self.create_object().ephemerides(quantities=17))

    def az_el(self):
        # The apparent azimuthal and elevation angle of the target.
        # Elevation is taken perpendicular to the local zentih direction.
        # UNITS: Degrees
        print(self.create_object().ephemerides(quantities=4))

    def R_ascdecl(self):
        # Apparent right ascension and declination of the
        # target with respect to the specified observing center/site.
        # UNITS: HMS/DMS or Degrees
        print(self.create_object().ephemerides(quantities=2))

    def angsep(self):
        # The angle between the center of a non-lunar target body and the
        # center of the primary body it revolves around, as seen by the observer.
        # UNITS: Arcseconds
        print(self.create_object().ephemerides(quantities=12))

    def target_angdia(self):
        # The angle subtended by the disk of the target seen by the observer, if it was fully illuminated.
        # UNITS: Arcseconds
        print(self.create_object().ephemerides(quantities=13))
