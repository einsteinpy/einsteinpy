import astropy.units as u

from einsteinpy import constant


@u.quantity_input(mass=u.kg)
def schwarzschild_radius(mass):
    """
    Schwarzschild radius is the radius defining the event horizon of a
    Schwarzschild black hole. It is characteristic radius associated with every
    quantity of mass.

    Parameters
    ----------
    mass : float
    """
    num = 2 * constant.G.value * mass
    deno = constant.c ** 2
    return num / deno

def scalar_factor(era, t, tuning_param=1):
    """
    Acceleration of the universe in cosmological models of Robertson Walker
    Flat Universe.

    Parameters
    ----------
    era : string
        Can be chosen from md (Matter Dominant),
        rd (Radiation Dominant) and ded (Dark Energy Dominant)
    t : float
        Time for the event (Value in range [0,1])
    tuning_param : float, optional
        Unit scaling factor, defaults to 1

    Returns
    -------
    a : float
        Value of scalar factor at time t.

    Raises
    ------
    ValueError : If era is not 'md' , 'rd', and 'ded'.

    """
    if era == 'md':
        return tuning_param * (t ** (2/3))
    elif era == 'rd':
        return tuning_param * (t ** (0.5))
    elif era == 'ded':
        hubble_const = (constant.Cosmo_Const / 3) ** 0.5
        val = np.e ** (hubble_const * t)
        return tuning_param * val
    else:
        raise ValueError('Passed era is out of syllabus!')
