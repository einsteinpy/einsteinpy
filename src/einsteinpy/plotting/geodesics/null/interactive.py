import os
import random

import numpy as np
from plotly import graph_objects as go
from plotly.offline import plot as saveplot


class InteractiveNullGeodesicPlotter:
    def __init__(self, bh_color="#ffcc00"):
        """
        Constructor

        Parameters
        ----------
        bh_color : string, optional
            Hexcpde for the color, used for the Black Hole
            Defaults to ``#ffcc00``

        """
        self.fig = go.Figure()
        self.bh_color = bh_color

        self._scene = dict(
            xaxis_title="X (GM/c^2)", yaxis_title="Y (GM/c^2)", zaxis_title="Z (GM/c^2)"
        )

    def _draw_bh(self, a):
        """
        Plots the Black Hole

        Parameters
        ----------
        a : float
            Black Hole Spin Parameter \
            0 <= a <= 1

        """
        r_outer = 1 + np.sqrt(1 - np.square(a))
        theta, phi = np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 50)
        THETA, PHI = np.meshgrid(theta, phi)

        X = r_outer * np.sin(PHI) * np.cos(THETA)
        Y = r_outer * np.sin(PHI) * np.sin(THETA)
        Z = r_outer * np.cos(PHI)

        self.fig.add_trace(
            go.Surface(
                x=X,
                y=Y,
                z=Z,
                name="Black Hole",
                opacity=0.75,
                showlegend=True,
                showscale=False,
            )
        )

    def plot(self, nullgeod, color="#{:06x}".format(random.randint(0, 0xFFFFFF))):
        """
        Plots the Null Geodesic

        Parameters
        ----------
        nullgeod : ~einsteinpy.geodesic.Nulllike
            Null Geodesic Object
        color : str, optional
            Hexcode for the color of the dashed lines, \
            that represent the Geodesic
            Picks a random color by default

        """
        a = nullgeod.a
        vals = nullgeod.trajectory[1]
        X = np.array([coord[1] for coord in vals])
        Y = np.array([coord[2] for coord in vals])
        Z = np.array([coord[3] for coord in vals])

        self._draw_bh(a)
        self.fig.add_trace(
            go.Scatter3d(
                x=X,
                y=Y,
                z=Z,
                name="Null Geodesic",
                marker=dict(size=1,),
                line=dict(color=color, width=2),
                showlegend=True,
            )
        )

    def show(self):
        """
        Shows plot during runtime

        Returns
        -------
        ~plotly.graph_objects.Figure

        """
        self.fig.update_layout(scene=self._scene)

        return self.fig

    def save(self, name="nullgeodesic.png"):
        """
        Saves plot locally

        Parameters
        ----------
        name : str, optional
            Name of the file, with extension

        """
        basename, ext = os.path.splitext(name)

        saveplot(self.fig, image=ext[1:], image_filename=basename, show_link=False)
