import random

from plotly.offline import init_notebook_mode

from .mplplotter import MatplotlibPlotter
from .plotlyplotter import PlotlyPlotter


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True
        if shell == "TerminalInteractiveShell":
            return False
        return False
    except NameError:
        return False


class Plotter(MatplotlibPlotter, PlotlyPlotter):
    """
    Class for automatically switching between Matplotlib and Plotly depending on platform used.
    """

    def __init__(self, attractor_color="#ffcc00"):
        """
        Constructor.

        Parameters
        ----------
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to #ffcc00.

        """
        if in_ipynb():
            PlotlyPlotter.__init__(self, attractor_color=attractor_color)
            self._notebook = True
            init_notebook_mode(connected=True)
        else:
            MatplotlibPlotter.__init__(self, attractor_color=attractor_color)
            self._notebook = False

    def plot(self, geodesic, color="#{:06x}".format(random.randint(0, 0xFFFFFF))):
        """

        Parameters
        ----------
        geodesic : ~einsteinpy.geodesic.Geodesic
            Geodesic of the body
        color : hex code RGB, optional
            Color of the dashed lines. Picks a random color by default.

        """
        if self._notebook:
            PlotlyPlotter.plot(self, geodesic, color=color)
        else:
            MatplotlibPlotter.plot(self, geodesic, color=color)

    def show(self):
        if self._notebook:
            return PlotlyPlotter.show(self)
        else:
            MatplotlibPlotter.show(self)

    def save(self, name):
        if self._notebook:
            PlotlyPlotter.save(self, name)
        else:
            MatplotlibPlotter.save(self, name)
