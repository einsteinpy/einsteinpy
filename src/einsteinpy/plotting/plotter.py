import random

from plotly.offline import init_notebook_mode

from .mplplotter import MatplotlibPlotter
from .plotlyplotter import PlotlyPlotter


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True
        elif shell == "TerminalInteractiveShell":
            return False
        else:
            return False
    except NameError:
        return False


class Plotter(MatplotlibPlotter, PlotlyPlotter):
    def __init__(self):
        if in_ipynb():
            PlotlyPlotter.__init__(self)
            self._notebook = True
            init_notebook_mode(connected=True)
        else:
            MatplotlibPlotter.__init__(self)
            self._notebook = False

    def plot(self, geodesic, color="#{:06x}".format(random.randint(0, 0xFFFFFF))):
        if self._notebook:
            PlotlyPlotter._draw_attractor(self)
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
