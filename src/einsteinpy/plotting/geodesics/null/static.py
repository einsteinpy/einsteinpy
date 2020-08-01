import random

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d


class StaticNullGeodesicPlotter:
    def __init__(self, ax=None, bh_color="#ffcc00", figsize=(10, 10)):
        """
        Constructor

        Parameters
        ----------
        ax: ~matplotlib.axes.Axes
            Matplotlib Axes object
        bh_color : string, optional
            Hexcpde for the color, used for the Black Hole
            Defaults to ``#ffcc00``
        figsize : tuple
            2-Tuple of Figure Size in inches
            Defaults to ``(10, 10)``

        """
        self.ax = ax
        if not self.ax:
            self.fig, self.ax = plt.subplots(figsize=figsize)
            self.ax = plt.axes(projection="3d")
            fontsize = max(figsize) + 3  # Scales well
            self.ax.set_xlabel("$X\\:(GM/c^2)$", fontsize=fontsize)
            self.ax.set_ylabel("$Y\\:(GM/c^2)$", fontsize=fontsize)
            self.ax.set_zlabel("$Z\\:(GM/c^2)$", fontsize=fontsize)

        self.bh_color = bh_color

    def _draw_bh(self, a, alpha=0.5):
        """
        Plots the Black Hole Outer Event Horizon

        Parameters
        ----------
        a : float
            Black Hole Spin Parameter
            0 <= a <= 1
        alpha : float, optional
            Opacity of the Black Hole Surface Plot
            0 <= alpha <= 1
            Defaults to ``0.5``

        """
        r_outer = 1 + np.sqrt(1 - np.square(a))
        theta, phi = np.linspace(0, 2 * np.pi, 50), np.linspace(0, np.pi, 50)
        THETA, PHI = np.meshgrid(theta, phi)

        X = r_outer * np.sin(PHI) * np.cos(THETA)
        Y = r_outer * np.sin(PHI) * np.sin(THETA)
        Z = r_outer * np.cos(PHI)

        self.ax.plot_surface(
            X,
            Y,
            Z,
            rstride=1,
            cstride=1,
            color=self.bh_color,
            antialiased=False,
            alpha=alpha,
            label="Black Hole",
        )

    def plot(self, nullgeod, color="#{:06x}".format(random.randint(0, 0xFFFFFF))):
        """
        Plots the Null Geodesic

        Parameters
        ----------
        nullgeod :  ~einsteinpy.geodesic.Nulllike
            Null Geodesic Object
        color : str, optional
            Hexcode for the color of the dashed lines, \
            that represent the Geodesic
            Picks a random color by default

        """
        a = nullgeod.a
        self._draw_bh(a)

        vals = nullgeod.trajectory[1]
        x = np.array(vals[:, 1])
        y = np.array(vals[:, 2])
        z = np.array(vals[:, 3])

        self.ax.plot(x, y, z, "--", color=color, label="Null Geodesic")

    def show(self, azim=-60, elev=30):
        """
        Adjusts the 3D view of the plot and \
        shows the plot during runtime

        Parameters
        ----------
        azim : float, optional
            Azimuthal viewing angle
            Defaults to ``-60`` Degrees

        elev : float, optional
            Elevation viewing angle
            Defaults to ``30`` Degrees

        """
        self.ax.view_init(azim=azim, elev=elev)
        plt.show()

    def save(self, name="nullgeodesic.png"):
        """
        Saves plot locally

        Parameters
        ----------
        name : str, optional
            Name of the file, with extension

        """
        plt.savefig(name)
