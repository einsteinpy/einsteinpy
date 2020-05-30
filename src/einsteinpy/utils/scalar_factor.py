import warnings

import astropy.units as u
import numpy as np

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
    warnings.warn("einteinpy.utils.scalar_factor.scalar_factor() \
        will be deprecated in Version 0.4. Please use \
        einsteinpy.metric.Metric.scalar_factor(t, era, tuning_param)!", \
        PendingDeprecationWarning) # - ????? Version/Warning
    
    T = t.to(u.s).value
    if era == "md":
        return tuning_param * (T ** (2 / 3))
    elif era == "rd":
        return tuning_param * (T ** (0.5))
    elif era == "ded":
        hubble_const = (constant.Cosmo_Const / 3) ** 0.5
        val = np.e ** (hubble_const.value * T)
        return tuning_param * val
    else:
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
    warnings.warn("einteinpy.utils.scalar_factor.scalar_factor_derivative() \
        will be deprecated in Version 0.4. Please use \
        einsteinpy.metric.Metric.scalar_factor_derivative(t, era, tuning_param)!", \
        PendingDeprecationWarning) # - ????? Version/Warning
    
    T = t.to(u.s).value
    if era == "md":
        return (2 / 3) * tuning_param * (T ** (-1 / 3))
    elif era == "rd":
        return 0.5 * tuning_param * (T ** (-0.5))
    elif era == "ded":
        hubble_const = (constant.Cosmo_Const / 3) ** 0.5
        val = hubble_const.value * (np.e ** (hubble_const.value * T))
        return tuning_param * val
    else:
        raise ValueError("Passed era should be either 'md', 'rd' or 'ded' ")
