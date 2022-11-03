import os
import random

import numpy as np
from plotly import graph_objects as go
from plotly.offline import plot as saveplot


class InteractiveGeodesicPlotter:
    def __init__(self, bh_colors=("#000", "#FFC"), draw_ergosphere=True):
        """
        Constructor

        Parameters
        ----------
        bh_colors : tuple, optional
            2-Tuple, containing hexcodes (Strings) for the colors,
            used for the Black Hole Event Horizon (Outer) and Ergosphere (Outer)
            Defaults to ``("#000", "#FFC")``
        draw_ergosphere : bool, optional
            Whether to draw the ergosphere
            Defaults to ``True``

        """
        self.fig = go.Figure()
        self.bh_colors = bh_colors
        self.draw_ergosphere = draw_ergosphere

    def _draw_bh(self, a):
        """
        Plots the Black Hole

        Parameters
        ----------
        a : float
            Dimensionless Spin Parameter of the Black Hole
            ``0 <= a <= 1``

        """
        colorscale_H = [(0, self.bh_colors[0]), (1, self.bh_colors[0])]
        colorscale_E = [(0, self.bh_colors[1]), (1, self.bh_colors[1])]

        theta, phi = np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 50)
        THETA, PHI = np.meshgrid(theta, phi)

        # Outer Event Horizon
        rh_outer = 1 + np.sqrt(1 - a**2)

        XH = rh_outer * np.sin(PHI) * np.cos(THETA)
        YH = rh_outer * np.sin(PHI) * np.sin(THETA)
        ZH = rh_outer * np.cos(PHI)

        self.fig.add_trace(
            go.Surface(
                x=XH,
                y=YH,
                z=ZH,
                name="BH Event Horizon (Outer)",
                opacity=0.75,
                showlegend=True,
                showscale=False,
                colorscale=colorscale_H,
            )
        )

        # Outer Ergosphere
        if self.draw_ergosphere:
            rE_outer = 1 + np.sqrt(1 - (a * np.cos(THETA) ** 2))

            XE = rE_outer * np.sin(PHI) * np.sin(THETA)
            YE = rE_outer * np.sin(PHI) * np.cos(THETA)
            ZE = rE_outer * np.cos(PHI)

            self.fig.add_trace(
                go.Surface(
                    x=XE,
                    y=YE,
                    z=ZE,
                    name="BH Ergosphere (Outer)",
                    opacity=0.5,
                    showlegend=True,
                    showscale=False,
                    colorscale=colorscale_E,
                )
            )

    def _draw_bh_2D(self, a, figsize=(6, 6)):
        """
        Plots the Black Hole in 2D

        Parameters
        ----------
        a : float
            Dimensionless Spin Parameter of the Black Hole
            ``0 <= a <= 1``
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(6, 6)``

        """
        theta = np.linspace(0, 2 * np.pi, 50)

        # Outer Event Horizon
        rh_outer = 1 + np.sqrt(1 - a**2)

        XH = rh_outer * np.sin(theta)
        YH = rh_outer * np.cos(theta)

        self.fig.add_trace(
            go.Scatter(
                x=XH,
                y=YH,
                fill="toself",
                name="BH Event Horizon (Outer)",
                opacity=0.75,
                showlegend=True,
                fillcolor=self.bh_colors[0],
            )
        )

        # Outer Ergosphere
        if self.draw_ergosphere:
            rE_outer = 1 + np.sqrt(1 - (a * np.cos(theta) ** 2))

            XE = rE_outer * np.sin(theta)
            YE = rE_outer * np.cos(theta)

            self.fig.add_trace(
                go.Scatter(
                    x=XE,
                    y=YE,
                    fill="toself",
                    name="BH Ergosphere (Outer)",
                    opacity=0.5,
                    showlegend=True,
                    fillcolor=self.bh_colors[1],
                )
            )

    def plot(
        self,
        geodesic,
        color="#{:06x}".format(random.randint(0, 0xFFFFFF)),
        title: str = "Geodesic Plot",
        aspect: str = "auto",
        aspect_ratio: dict = dict(x=1, y=1, z=1),
    ):
        """
        Plots the Geodesic

        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        color : str, optional
            Hexcode (String) for the color of the
            dashed lines, that represent the Geodesic
            Picks a random color by default
        title : str, optional
            Plot title
        aspect : {"auto", "data", "manual", "cube"}
            Aspect ratio for plot axes
            Defaults to "auto"
        aspect_ratio : dict, optional
            Aspect ratio to define if aspect type is "manual"
            Defaults to ``dict(x=1, y=1, z=1)``

        Raises
        ------
        ValueError
            If ``aspect`` does not take values from ``{"auto", "data", "manual", "cube"}``

        """
        aspects = ["auto", "data", "manual", "cube"]

        if aspect not in aspects:
            raise ValueError(
                f"Invalid aspect type. Expected one of {aspects}. Received '{aspect}'."
            )

        traj = geodesic.trajectory[1]
        X = traj[:, 1]
        Y = traj[:, 2]
        Z = traj[:, 3]

        a = geodesic.metric_params[0]
        self._draw_bh(a)
        self.fig.add_trace(
            go.Scatter3d(
                x=X,
                y=Y,
                z=Z,
                name=geodesic.kind + " Geodesic",
                marker=dict(
                    size=1,
                ),
                line=dict(color=color, width=2),
                showlegend=True,
            )
        )

        _title_3D = dict(text=title, y=0.9, x=0.46, xanchor="center", yanchor="top")

        _scene_3D = dict(
            xaxis_title="X (GM/c^2)", yaxis_title="Y (GM/c^2)", zaxis_title="Z (GM/c^2)"
        )

        self.fig.update_layout(
            title=_title_3D,
            scene=_scene_3D,
            scene_aspectmode=aspect,
            scene_aspectratio=aspect_ratio,
        )

    def plot2D(
        self,
        geodesic,
        coordinates=(1, 2),
        figsize=(6, 6),
        color="#{:06x}".format(random.randint(0, 0xFFFFFF)),
        title: str = "",
    ):
        """
        Plots the Geodesic in 2D

        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        coordinates : tuple, optional
            2-Tuple, containing labels for coordinates & momenta to plot
            Labels for ``X1, X2, X3, P1, P2, P3`` are ``(1, 2, 3, 4, 5, 6)``
            Defaults to ``(1, 2)`` (X, Y)
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(6, 6)``
        color : str, optional
            Hexcode (String) for the color of the
            dashed lines, that represent the Geodesic
            Picks a random color by default
        title : str, optional
            Plot title

        Raises
        ------
        IndexError
            If indices in ``coordinates`` do not take values from ``(1, 2, 3)``

        """
        a = geodesic.metric_params[0]
        self._draw_bh_2D(a, figsize)

        traj = geodesic.trajectory[1]
        A = coordinates[0]
        B = coordinates[1]

        if A not in (1, 2, 3) or B not in (1, 2, 3):
            raise IndexError(
                """
                Please ensure that indices in `coordinates` take two of these values: `(1, 2, 3)`.
                E.g., Indices for `X1, X2, X3` are `(1, 2, 3)`.
                """
            )

        self.fig.add_trace(
            go.Scatter(
                x=traj[:, A],
                y=traj[:, B],
                line=dict(color=color, width=2),
                name=geodesic.kind + " Geodesic",
                mode="lines",
            )
        )

        self.fig.update_layout(
            xaxis_title=f"X{coordinates[0]} (GM/c^2)",
            yaxis_title=f"X{coordinates[1]} (GM/c^2)",
        )

        if title:
            self.fig.layout.title = {
                "text": title,
                "x": 0.46,
                "xanchor": "center",
                "y": 0.9,
                "yanchor": "top",
            }

    def parametric_plot(
        self,
        geodesic,
        colors=("#00FFFF", "#FF00FF", "#FFFF00"),
        title: str = "Parametric Plot",
    ):
        """
        Plots the coordinates of the Geodesic, against Affine Parameter

        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        colors : tuple, optional
            3-Tuple, containing hexcodes (Strings) for the color
            of the lines, for each of the 3 coordinates
            Defaults to ``("#00FFFF", "#FF00FF", "#FFFF00")``
        title : str, optional
            Plot title

        """
        coords = geodesic.coords
        traj = geodesic.trajectory
        lambdas = traj[0]
        X1 = traj[1][:, 1]
        X2 = traj[1][:, 2]
        X3 = traj[1][:, 3]

        self.fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=X1,
                name=f"X1 ({coords})",
                line=dict(color=colors[0], width=2),
                mode="lines",
            )
        )
        self.fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=X2,
                name=f"X2 ({coords})",
                line=dict(color=colors[1], width=2),
                mode="lines",
            )
        )
        self.fig.add_trace(
            go.Scatter(
                x=lambdas,
                y=X3,
                name=f"X3 ({coords})",
                line=dict(color=colors[2], width=2),
                mode="lines",
            )
        )

        _title_2D = dict(text=title, y=0.9, x=0.46, xanchor="center", yanchor="top")

        self.fig.update_layout(
            title=_title_2D,
            xaxis_title="Lambda (Affine Parameter)",
            yaxis_title="Coordinates",
        )

    def show(self):
        """
        Shows plot during runtime

        Returns
        -------
        ~plotly.graph_objects.Figure

        """
        return self.fig

    def clear(self):
        """
        Clears plot during runtime

        """
        self.fig.data = []
        self.fig.layout = {}

    def save(self, name="Geodesic.png"):
        """
        Saves plot locally

        Parameters
        ----------
        name : str, optional
            Name of the file, with extension
            Defaults to ``Geodesic.png``

        """
        basename, ext = os.path.splitext(name)

        saveplot(self.fig, image=ext[1:], image_filename=basename, show_link=False)
