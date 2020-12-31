from astropy import units as u

from einsteinpy.constant import G, c

__all__ = ["primitive"]


def primitive(*args):
    """
    Strips out units and returns numpy.float64 values \
    out of ``astropy.units.quantity.Quantity``

    Parameters
    ----------
    *args : iterable
        ``astropy.units.quantity.Quantity`` objects, who ``value`` is required

    Returns
    -------
    primitive_args : list
        List of ``numpy.float64`` values, obtained from ``Quantity`` objects

    """
    primitive_args = []
    for item in args:
        if isinstance(item, u.Quantity):
            primitive_args.append(item.value)
        else:
            primitive_args.append(item)

    return primitive_args
