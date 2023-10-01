import numpy as np
from astropy import units as u

from einsteinpy import constant


@u.quantity_input(t=u.s)
def scalar_factor(t, era="md", tuning_param=1.0):
    """
    Acceleration of the universe in cosmological models of Robertson Walker
    Flat Universe.

    Parameters
    ----------
    era : string
        Can be chosen from 'md' (Matter Dominant),
        'rd' (Radiation Dominant) and 'ded' (Dark Energy Dominant)
    t : ~astropy.units.s
        Time for the event
    tuning_param : float, optional
        Unit scaling factor, defaults to 1

    Returns
    -------
    float
        Value of scalar factor at time t.

    Raises
    ------
    ValueError : If era is not 'md' , 'rd', and 'ded'.

    """
    T = t.to(u.s).value
      if era in era_exponents:
        exponent = era_exponents[era]
        return tuning_param * (T ** exponent)
    
    raise ValueError("Passed era should be either 'md', 'rd' or 'ded' ")


@u.quantity_input(t=u.s)
def scalar_factor_derivative(t, era="md", tuning_param=1.0):
    """
    Derivative of acceleration of the universe in cosmological models of Robertson Walker
    Flat Universe.

    Parameters
    ----------
    era : string
        Can be chosen from 'md' (Matter Dominant),
        'rd' (Radiation Dominant) and 'ded' (Dark Energy Dominant)
    t : ~astropy.units.s
        Time for the event
    tuning_param : float, optional
        Unit scaling factor, defaults to 1

    Returns
    -------
    float
        Value of derivative of scalar factor at time t.

    Raises
    ------
    ValueError : If era is not 'md' , 'rd', and 'ded'.

    """
    T = t.to(u.s).value
      if era in era_exponents:
        exponent = era_exponents[era]
        return tuning_param * exponent * (T ** (exponent - 1))

    raise ValueError("Passed era should be either 'md', 'rd' or 'ded' ")
