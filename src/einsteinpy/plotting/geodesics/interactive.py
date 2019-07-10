import os
import random

import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot as saveplot


class InteractiveGeodesicPlotter:
    def __init__(self, attractor_color="#ffcc00"):
        """
        Constructor.

        Parameters
        ----------
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to #ffcc00.

        """
        self.fig = go.FigureWidget()
        self.attractor_color = attractor_color
        self.attractor_present = False
        self._layout = go.Layout(
            autosize=True,
            xaxis=dict(title="x (m)", constrain="domain"),
            yaxis=dict(title="y (m)", scaleanchor="x"),
        )

    def _draw_attractor(self, radius, xarr, yarr):
        self.attractor_present = True
        self.fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="markers",
                name="attractor",
                marker=dict(size=10, color=self.attractor_color, line=dict(width=2)),
            )
        )

    def plot(self, geodesic, color="#{:06x}".format(random.randint(0, 0xFFFFFF))):
        """

        Parameters
        ----------
        geodesic : ~einsteinpy.geodesic.Geodesic
            Geodesic of the body
        color : hex code RGB, optional
            Color of the dashed lines. Picks a random color by default.

        """
        vals = geodesic.trajectory
        x = np.array([coord[1] for coord in vals])
        y = np.array([coord[2] for coord in vals])

        if not self.attractor_present:
            self._draw_attractor(geodesic.metric.scr, x, y)
        self.fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                name="geodesic",
                marker=dict(size=5, color=color, line=dict(width=2)),
            )
        )

    def show(self):
        """
        Method to show plots during runtime.

        Returns
        -------
        ~plotly.graph_objs.FigureWidget
        """
        self.fig.layout.update(self._layout)
        return self.fig

    def save(self, name="geodesic.png"):
        """
        Method to save plots locally.

        Parameters
        ----------
        name : str, optional
            Name of the file with extension.
        """
        basename, ext = os.path.splitext(name)
        saveplot(self.fig, image=ext[1:], image_filename=basename, show_link=False)
