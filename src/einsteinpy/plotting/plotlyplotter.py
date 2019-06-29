import os
import random

import plotly.graph_objs as go
import numpy as np
from plotly.offline import plot as saveplot


class PlotlyPlotter:
    def __init__(self, attractor_color="#ffcc00"):
        self.fig = go.FigureWidget()
        self.attractor_color = attractor_color
        self.attractor_present = False

    def _draw_attractor(self):
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
        vals = geodesic.trajectory
        r = np.array([coord[1] for coord in vals])
        phi = np.array([coord[3] for coord in vals])

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        if not self.attractor_present:
            self._draw_attractor()
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
        return self.fig

    def save(self, name="geodesic.png"):
        basename, ext = os.path.splitext(name)
        saveplot(
            self.fig,
            image=ext[1:],
            image_filename=basename,
            show_link=False,
        )
