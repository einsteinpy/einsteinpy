import random
import warnings

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d


class StaticGeodesicPlotter:
    def __init__(self, ax=None, bh_colors=("#000", "#FFC"), draw_ergosphere=True):
        """
        Constructor

        Parameters
        ----------
        ax: ~matplotlib.axes.Axes
            Matplotlib Axes object
            To be deprecated in Version 0.5.0
            Since Version 0.4.0, `StaticGeodesicPlotter`
            automatically creates a new Axes Object.
            Defaults to ``None``
        bh_colors : tuple, optional
            2-Tuple, containing hexcodes (Strings) for the colors,
            used for the Black Hole Event Horizon (Outer) and Ergosphere (Outer)
            Defaults to ``("#000", "#FFC")``
        draw_ergosphere : bool, optional
            Whether to draw the ergosphere
            Defaults to `True`

        """
        self.ax = ax
        self.bh_colors = bh_colors
        self.draw_ergosphere = draw_ergosphere

        if ax is not None:
            warnings.warn(
                """
                Argument `ax` will be removed in Version 0.5.0.
                Since Version 0.4.0, `StaticGeodesicPlotter` automatically
                creates a new Axes Object.
                """,
                PendingDeprecationWarning,
            )

    def _draw_bh(self, a, figsize=(6, 6)):
        """
        Plots the Black Hole in 3D

        Parameters
        ----------
        a : float
            Dimensionless Spin Parameter of the Black Hole
            ``0 <= a <= 1``
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(6, 6)``

        """
        self.fig, self.ax = plt.subplots(figsize=figsize)
        fontsize = max(figsize) + 3
        self.fig.set_size_inches(figsize)
        self.ax = plt.axes(projection="3d")
        self.ax.set_xlabel("$X\\:(GM/c^2)$", fontsize=fontsize)
        self.ax.set_ylabel("$Y\\:(GM/c^2)$", fontsize=fontsize)
        self.ax.set_zlabel("$Z\\:(GM/c^2)$", fontsize=fontsize)

        theta, phi = np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 50)
        THETA, PHI = np.meshgrid(theta, phi)

        # Outer Event Horizon
        rh_outer = 1 + np.sqrt(1 - a**2)

        XH = rh_outer * np.sin(PHI) * np.cos(THETA)
        YH = rh_outer * np.sin(PHI) * np.sin(THETA)
        ZH = rh_outer * np.cos(PHI)

        surface1 = self.ax.plot_surface(
            XH,
            YH,
            ZH,
            rstride=1,
            cstride=1,
            color=self.bh_colors[0],
            antialiased=False,
            alpha=0.2,
            label="BH Event Horizon (Outer)",
        )

        surface1._facecolors2d = surface1._facecolor3d
        surface1._edgecolors2d = surface1._edgecolor3d

        # Outer Ergosphere
        if self.draw_ergosphere:
            rE_outer = 1 + np.sqrt(1 - (a * np.cos(THETA) ** 2))

            XE = rE_outer * np.sin(PHI) * np.sin(THETA)
            YE = rE_outer * np.sin(PHI) * np.cos(THETA)
            ZE = rE_outer * np.cos(PHI)

            surface2 = self.ax.plot_surface(
                XE,
                YE,
                ZE,
                rstride=1,
                cstride=1,
                color=self.bh_colors[1],
                antialiased=False,
                alpha=0.1,
                label="BH Ergosphere (Outer)",
            )

            surface2._facecolors2d = surface2._facecolor3d
            surface2._edgecolors2d = surface2._edgecolor3d

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
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.fig.set_size_inches(figsize)

        theta = np.linspace(0, 2 * np.pi, 50)

        # Outer Event Horizon
        rh_outer = 1 + np.sqrt(1 - a**2)

        XH = rh_outer * np.sin(theta)
        YH = rh_outer * np.cos(theta)

        self.ax.fill(
            XH, YH, self.bh_colors[0], alpha=0.2, label="BH Event Horizon (Outer)"
        )

        # Outer Ergosphere
        if self.draw_ergosphere:
            rE_outer = 1 + np.sqrt(1 - (a * np.cos(theta) ** 2))

            XE = rE_outer * np.sin(theta)
            YE = rE_outer * np.cos(theta)

            self.ax.fill(
                XE, YE, self.bh_colors[1], alpha=0.1, label="BH Ergosphere (Outer)"
            )

    def plot(
        self,
        geodesic,
        figsize=(6, 6),
        color="#{:06x}".format(random.randint(0, 0xFFFFFF)),
        title: str = "",
        aspect: str = "auto",
    ):
        """
        Plots the Geodesic

        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(6, 6)``
        color : str, optional
            Hexcode (String) for the color of the
            dashed lines, that represent the Geodesic
            Picks a random color by default
        title : str, optional
            Plot title
        aspect : {"auto", "equal", "equalxy", "equalyz", "equalxz"}
            Aspect ratio for plot axes
            Defaults to "auto"

        Raises
        ------
        ValueError
            If ``aspect`` does not take values from ``{"auto", "equal", "equalxy", "equalyz", "equalxz"}``

        """
        aspects = ["auto", "equal", "equalxy", "equalyz", "equalxz"]

        if aspect not in aspects:
            raise ValueError(
                f"Invalid aspect type. Expected one of {aspects}. Received '{aspect}'."
            )

        a = geodesic.metric_params[0]
        self._draw_bh(a, figsize)

        traj = geodesic.trajectory[1]
        x = traj[:, 1]
        y = traj[:, 2]
        z = traj[:, 3]

        self.ax.plot(x, y, z, "--", color=color, label=geodesic.kind + " Geodesic")
        self.ax.set_aspect(aspect)

        if title:
            self.ax.set_title(title)

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
            2-Tuple, containing labels for coordinates to plot
            Labels for ``X1, X2, X3`` are ``(1, 2, 3)``
            Defaults to ``(1, 2)`` (X, Y)
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(6, 6)``
        color : str, optional
            Hexcode (String) for the color of the
            dashed lines, that represent the Geodesic
            Picks a random color by default
        title: str, optional
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
                Please ensure, that indices in `coordinates` take two of these values: `(1, 2, 3)`.
                Indices for `X1, X2, X3` are `(1, 2, 3)`.
                """
            )

        fontsize = max(figsize) + 3
        self.ax.set_xlabel(f"$X{coordinates[0]}\\:(GM/c^2)$", fontsize=fontsize)
        self.ax.set_ylabel(f"$X{coordinates[1]}\\:(GM/c^2)$", fontsize=fontsize)

        self.ax.plot(
            traj[:, A], traj[:, B], "--", color=color, label=geodesic.kind + " Geodesic"
        )

        if title:
            self.ax.set_title(title)

    def parametric_plot(
        self,
        geodesic,
        figsize=(8, 6),
        colors=("#00FFFF", "#FF00FF", "#FFFF00"),
        title: str = "",
    ):
        """
        Plots the coordinates of the Geodesic, against Affine Parameter

        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        figsize : tuple, optional
            2-Tuple of Figure Size in inches
            Defaults to ``(8, 6)``
        colors : tuple, optional
            3-Tuple, containing hexcodes (Strings) for the color
            of the lines, for each of the 3 coordinates
            Defaults to ``("#00FFFF", "#FF00FF", "#00FFFF")``
        title : str, optional
            Plot title

        """
        self.fig, self.ax = plt.subplots(figsize=figsize)
        fontsize = max(figsize) + 3
        self.fig.set_size_inches(figsize)
        self.ax.set_xlabel(r"Affine Paramter, $\lambda$", fontsize=fontsize)
        self.ax.set_ylabel("Coordinates", fontsize=fontsize)

        coords = geodesic.coords
        traj = geodesic.trajectory
        lambdas = traj[0]
        X1 = traj[1][:, 1]
        X2 = traj[1][:, 2]
        X3 = traj[1][:, 3]

        self.ax.plot(lambdas, X1, color=colors[0], label=f"X1 ({coords})")
        self.ax.plot(lambdas, X2, color=colors[1], label=f"X2 ({coords})")
        self.ax.plot(lambdas, X3, color=colors[2], label=f"X3 ({coords})")

        self.ax.set_title(title)

    def animate(
        self, geodesic, interval=10, color="#{:06x}".format(random.randint(0, 0xFFFFFF))
    ):
        """
        Parameters
        ----------
        geodesic : einsteinpy.geodesic.*
            Geodesic Object
        interval : int, optional
            Time (in milliseconds) between frames
            Defaults to ``10``
        color : str, optional
            Hexcode (String) for the color of the
            dashed lines, that represent the Geodesic
            Picks a random color by default

        """
        a = geodesic.metric_params[0]
        self._draw_bh(a)

        traj = geodesic.trajectory
        x = traj[1][:, 1]
        y = traj[1][:, 2]
        z = traj[1][:, 3]
        N = x.shape[0]

        x_max, x_min = max(x), min(x)
        y_max, y_min = max(y), min(y)
        z_max, z_min = max(z), min(z)
        margin_x = (x_max - x_min) * 0.2
        margin_y = (y_max - y_min) * 0.2
        margin_z = (z_max - z_min) * 0.2

        self.ax.set_xlim3d([x_min - margin_x, x_max + margin_x])
        self.ax.set_ylim3d([y_min - margin_y, y_max + margin_y])
        self.ax.set_zlim3d([z_min - margin_z, z_max + margin_z])

        data = traj[1][:, 1:4].T
        (line,) = self.ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

        def _update(num, data, line):
            line.set_data(data[:2, :num])
            line.set_3d_properties(data[2, :num])

            return (line,)

        self.ani = FuncAnimation(
            self.fig, _update, N, fargs=(data, line), interval=interval, blit=True
        )

    def show(self, azim=-60, elev=30):
        """
        Adjusts the 3D view of the plot and \
        shows the plot during runtime. For Parametric Plots,
        only the plot is displayed.

        Parameters
        ----------
        azim : float, optional
            Azimuthal viewing angle
            Defaults to ``-60`` Degrees

        elev : float, optional
            Elevation viewing angle
            Defaults to ``30`` Degrees

        """
        figsize = self.fig.get_size_inches()
        fontsize = max(figsize) + 1.5
        if self.ax.name == "3d":
            self.ax.view_init(azim=azim, elev=elev)
        plt.legend(prop={"size": fontsize})

        plt.show()

    def clear(self):
        """
        Clears plot during runtime

        """
        self.fig.clf()

    def save(self, name="Geodesic.png"):
        """
        Saves plot locally
        Should be called before ``show()``, as
        ``show()`` erases current figure's contents.

        Parameters
        ----------
        name : str, optional
            Name of the file, with extension
            Defaults to ``Geodesic.png``

        """
        if self.ax.name != "3d" and name == "Geodesic.png":
            name = "Parametric.png"

        plt.savefig(name)
