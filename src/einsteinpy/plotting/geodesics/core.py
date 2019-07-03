import random

from .interactive import InteractiveGeodesicPlotter
from .static import StaticGeodesicPlotter


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__  # type: ignore
        if shell == "ZMQInteractiveShell":
            return InteractiveGeodesicPlotter
        if shell == "TerminalInteractiveShell":
            return StaticGeodesicPlotter
        return StaticGeodesicPlotter
    except NameError:
        return StaticGeodesicPlotter


class GeodesicPlotter(in_ipynb()):  # type: ignore
    """
    Class for automatically switching between Matplotlib and Plotly depending on platform used.
    """
