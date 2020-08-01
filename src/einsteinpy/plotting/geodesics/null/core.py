import random

from .interactive import InteractiveNullGeodesicPlotter
from .static import StaticNullGeodesicPlotter


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__  # type: ignore
        if shell == "ZMQInteractiveShell":
            return InteractiveNullGeodesicPlotter
        if shell == "TerminalInteractiveShell":
            return StaticNullGeodesicPlotter
        return StaticNullGeodesicPlotter
    except NameError:
        return StaticNullGeodesicPlotter


class NullGeodesicPlotter(in_ipynb()):  # type: ignore
    """
    Class for automatically switching between Matplotlib and Plotly depending on platform used.
    """
