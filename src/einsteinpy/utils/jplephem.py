from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric

"""
The velocity cannot be calculated for the Moon.
To just get the position, use get_body_helio()
"""


def _get_body_helio_posvel(body, time, ephemeris=None, get_velocity=True):
    """
    Function to get velocity and position of a solar system body w.r.t the sun.

    Parametes
    ---------
    body: str or other
        The solar system body for which to calculate positions.
        Can also be a kernel specifier (list of 2-tuples) if the ephemeris is a JPL kernel.

    time:
    Time of observation.

    ephemeris : str, optional
    Ephemeris to use.
    By default, use the one set with astropy.coordinates.solar_system_ephemeris.set
    """
    if get_velocity:
        body_pos_bary, body_vel_bary = get_body_barycentric_posvel(
            body, time, ephemeris
        )
        sun_pos_bary, sun_vel_bary = get_body_barycentric_posvel(
            "sun", time, ephemeris
        )

        body_pos_helio = body_pos_bary - sun_pos_bary
        body_vel_helio = body_vel_bary - sun_vel_bary
    else:
        body_pos_bary = get_body_barycentric(body, time, ephemeris)
        sun_pos_bary = get_body_barycentric("sun", time, ephemeris)

        body_pos_helio = body_pos_bary - sun_pos_bary

    return (body_pos_helio, body_vel_helio) if get_velocity else body_pos_helio


def get_body_helio_posvel(body, time, ephemeris=None):
    return _get_body_helio_posvel(body, time, ephemeris)


def get_body_helio(body, time, ephemeris=None):
    return _get_body_helio_posvel(body, time, ephemeris, get_velocity=False)
