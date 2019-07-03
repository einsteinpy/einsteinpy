import random

from plotly.offline import init_notebook_mode

from .mplplotter import MatplotlibPlotter
from .plotlyplotter import PlotlyPlotter


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__  # type: ignore
        if shell == "ZMQInteractiveShell":
            return PlotlyPlotter
        if shell == "TerminalInteractiveShell":
            return MatplotlibPlotter
        return MatplotlibPlotter
    except NameError:
        return MatplotlibPlotter


class Plotter(in_ipynb()):  # type: ignore
    """
    Class for automatically switching between Matplotlib and Plotly depending on platform used.
    """
