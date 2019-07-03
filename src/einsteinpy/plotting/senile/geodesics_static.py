import random

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from einsteinpy.utils import schwarzschild_radius


class StaticGeodesicPlotter:
    """
    Class for plotting static matplotlib plots and animations.
    """

    def __init__(
        self,
        time=0 * u.s,
        ax=None,
        attractor_radius_scale=-1.0,
        attractor_color="#ffcc00",
    ):
        """
        Constructor.

        Parameters
        ----------
        time : ~astropy.units.quantity.Quantity
            Time of start, defaults to 0 seconds.
        attractor_radius_scale : float, optional
            Scales the attractor radius by the value given. Default is 1. It is used to make plots look more clear if needed.
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to #ffcc00.

        """
        self.ax = ax
        if not self.ax:
            self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.time = time
        self._attractor_present = False
        self.attractor_color = attractor_color
        self.attractor_radius_scale = attractor_radius_scale
        self.__xarr = np.array([])
        self.__yarr = np.array([])
        self.get_curr_plot_radius = 0

    def plot_trajectory(self, geodesic, color, only_points=False):
        """

        Parameters
        ----------
        geodesic : ~einsteinpy.geodesic.Geodesic
            Geodesic of the body.
        color : string
            Color of the Geodesic

        """
        vals = geodesic.trajectory

        # time = np.array([coord[0] for coord in vals])
        r = np.array([coord[1] for coord in vals])
        # theta = np.array([coord[2] for coord in vals])
        phi = np.array([coord[3] for coord in vals])

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        self.__xarr = x
        self.__yarr = y

        if only_points:
            return x, y

        lines = self.ax.plot(x, y, "--", color=color)

        return lines, x[-1], y[-1]

    def plot_attractor(self):
        if not self._attractor_present:
            self._draw_attractor()

    def mindist(self, x, y):
        return np.sqrt(x * x + y * y)

    def _draw_attractor(self):
        radius = schwarzschild_radius(self.mass)
        if self.attractor_radius_scale == -1.0:
            minrad_nooverlap = self.mindist(self.__xarr[0], self.__yarr[0])
            for i in range(0, len(self.__xarr)):
                minrad_nooverlap = min(
                    minrad_nooverlap, self.mindist(self.__xarr[i], self.__yarr[i])
                )

            xlen = max(self.__xarr) - min(self.__xarr)
            ylen = max(self.__yarr) - min(self.__yarr)
            minlen_plot = min(xlen, ylen)
            mulitplier = minlen_plot / (12 * radius)
            min_radius = radius * mulitplier

            radius = min(min_radius, minrad_nooverlap)
            self.get_curr_plot_radius = radius
            self.ax.add_patch(
                mpl.patches.Circle((0, 0), radius, lw=0, color=self.attractor_color)
            )
        else:
            radius = radius * self.attractor_radius_scale
            self.get_curr_plot_radius = radius
            self.ax.add_patch(
                mpl.patches.Circle(
                    (0, 0), radius.value, lw=0, color=self.attractor_color
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

        Returns
        -------
        lines : list
            A list of Line2D objects representing the plotted data.

        """
        _, _ = self.plot_trajectory(geodesic, only_points=True, color=color)
        self.mass = geodesic.attractor.mass
        self.plot_attractor()
        self._attractor_present = True

        lines, x0, y0 = self.plot_trajectory(geodesic, color)

        l, = self.ax.plot(x0, y0, "o", mew=0, color=lines[0].get_color())
        lines.append(l)

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines

    def animate(
        self, geodesic, color="#{:06x}".format(random.randint(0, 0xFFFFFF)), interval=50
    ):
        """

        Parameters
        ----------
        geodesic : ~einsteinpy.geodesic.Geodesic
            Geodesic of the body.
        color : hex code RGB, optional
            Color of the dashed lines. Picks a random color by default.
        interval : int, optional
            Control the time between frames. Add time in milliseconds.

        """

        pos_x, pos_y = self.plot_trajectory(geodesic, only_points=True, color=color)
        self.mass = geodesic.attractor.mass
        x_max, x_min = max(pos_x), min(pos_x)
        y_max, y_min = max(pos_y), min(pos_y)
        margin_x = (x_max - x_min) * 0.1
        margin_y = (y_max - y_min) * 0.2
        frames = pos_x.shape[0]

        pic, = self.ax.plot([], [], "--", color=color)

        plt.xlim(x_min - margin_x, x_max + margin_x)
        plt.ylim(y_min - margin_y, y_max + margin_y)

        self.plot_attractor()
        self._attractor_present = True

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        def _update(frame):
            pic.set_xdata(pos_x[: frame + 1])
            pic.set_ydata(pos_y[: frame + 1])
            return (pic,)

        self.ani = FuncAnimation(
            self.fig, _update, frames=frames, interval=interval, blit=True
        )

    def save(self, name="static_geodesic.png"):
        plt.savefig(name)

    def show(self):
        plt.show()
