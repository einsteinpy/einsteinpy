import os
import random

import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot as saveplot


class InteractiveGeodesicPlotter:
    def __init__(
        self, attractor_radius_scale=-1.0, attractor_color="#ffcc00", use_3d=False
    ):
        """
        Constructor.

        Parameters
        ----------
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to #ffcc00.

        """
        self.fig = go.Figure()
        self.attractor_radius_scale = attractor_radius_scale
        self.use_3d = use_3d
        self.attractor_color = attractor_color
        self.attractor_present = False
        if use_3d:
            self._layout = go.Layout(
                autosize=True,
                scene=dict(
                    xaxis=dict(title="x (m)"),
                    yaxis=dict(title="y (m)"),
                    zaxis=dict(title="z (m)"),
                ),
            )
        else:
            self._layout = go.Layout(
                autosize=True,
                xaxis=dict(title="x (m)", constrain="domain"),
                yaxis=dict(title="y (m)", scaleanchor="x"),
            )

    def _mindist(self, x, y, z=0):
        return np.sqrt(x * x + y * y + z * z)

    def _draw_attractor(self, radius, xarr, yarr):
        self.attractor_present = True
        if self.attractor_radius_scale == -1.0:
            minrad_nooverlap = self._mindist(xarr[0], yarr[0])
            for i, _ in enumerate(xarr):
                minrad_nooverlap = min(
                    minrad_nooverlap, self._mindist(xarr[i], yarr[i])
                )

            xlen = max(xarr) - min(xarr)
            ylen = max(yarr) - min(yarr)
            minlen_plot = min(xlen, ylen)
            min_radius = minlen_plot / 12
            radius = min(min_radius, minrad_nooverlap)
        else:
            radius = radius.value * self.attractor_radius_scale

        if self.use_3d:
            self.fig.add_trace(
                go.Scatter3d(
                    x=[0],
                    y=[0],
                    z=[0],
                    mode="markers",
                    name="attractor",
                    marker=dict(size=0, color=self.attractor_color, line=dict(width=0)),
                )
            )
        else:
            self.fig.add_trace(
                go.Scatter(
                    x=[0],
                    y=[0],
                    mode="markers",
                    name="attractor",
                    marker=dict(size=0, color=self.attractor_color),
                )
            )
            self._layout.shapes = [
                go.layout.Shape(
                    type="circle",
                    xref="x",
                    yref="y",
                    fillcolor=self.attractor_color,
                    x0=-radius,
                    y0=-radius,
                    x1=radius,
                    y1=radius,
                    line_width=0,
                )
            ]

    def _set_scaling(self, x_range, y_range, z_range, lim):
        if x_range < lim and y_range < lim and z_range < lim:
            return
        if x_range < lim:
            self._layout.scene.xaxis.range = [-lim, lim]
        if y_range < lim:
            self._layout.scene.yaxis.range = [-lim, lim]
        if z_range < lim:
            self._layout.scene.zaxis.range = [-lim, lim]

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
        x = np.array(vals[:, 1])
        y = np.array(vals[:, 2])

        if not self.attractor_present:
            self._draw_attractor(geodesic.metric.scr, x, y)
        if self.use_3d:
            z = np.array(vals[:, 3])
            x_range = max(x) - min(x)
            y_range = max(y) - min(y)
            z_range = max(z) - min(z)
            self._set_scaling(x_range, y_range, z_range, 1e-5)
            self.fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="lines",
                    name="geodesic",
                    marker=dict(size=5, color=color, line=dict(width=2)),
                )
            )
        else:
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
        ~plotly.graph_objects.Figure
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
