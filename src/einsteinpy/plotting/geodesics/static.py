import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle


class StaticGeodesicPlotter:
    def __init__(self, ax=None, attractor_radius_scale=-1.0, attractor_color="#ffcc00"):
        """
        Constructor.

        Parameters
        ----------
        ax: ~matplotlib.axes.Axes
            Matplotlib Axes object.
        attractor_radius_scale : float, optional
            Scales the attractor radius by the value given. Default is 1.
            It is used to make plots look more clear if needed.
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to #ffcc00.

        """
        self.ax = ax
        if not self.ax:
            self.fig, self.ax = plt.subplots(figsize=(6, 6))
            self.ax.set_xlabel("$x$ (km)")
            self.ax.set_ylabel("$y$ (km)")
            self.ax.set_aspect(1)
        self.attractor_radius_scale = attractor_radius_scale
        self.attractor_color = attractor_color
        self.attractor_present = False

    def _mindist(self, x, y):
        return np.sqrt(x * x + y * y)

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
            mulitplier = minlen_plot / (12 * radius)
            min_radius = radius * mulitplier

            radius = min(min_radius, minrad_nooverlap)
            self.get_curr_plot_radius = radius
            self.ax.add_patch(Circle((0, 0), radius, lw=0, color=self.attractor_color))
        else:
            radius = radius * self.attractor_radius_scale
            self.get_curr_plot_radius = radius
            self.ax.add_patch(
                Circle((0, 0), radius.value, lw=0, color=self.attractor_color)
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
        self.ax.plot(x, y, "--", color=color)
        self.ax.plot(x[-1], y[-1], "o", mew=0, color=color)

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

        vals = geodesic.trajectory
        x = np.array([coord[1] for coord in vals])
        y = np.array([coord[2] for coord in vals])
        x_max, x_min = max(x), min(x)
        y_max, y_min = max(y), min(y)
        margin_x = (x_max - x_min) * 0.1
        margin_y = (y_max - y_min) * 0.2
        frames = x.shape[0]

        pic, = self.ax.plot([], [], "--", color=color)

        plt.xlim(x_min - margin_x, x_max + margin_x)
        plt.ylim(y_min - margin_y, y_max + margin_y)

        if not self.attractor_present:
            self._draw_attractor(geodesic.metric.scr, x, y)

        def _update(frame):
            pic.set_xdata(x[: frame + 1])
            pic.set_ydata(y[: frame + 1])
            return (pic,)

        self.ani = FuncAnimation(
            self.fig, _update, frames=frames, interval=interval, blit=True
        )

    def show(self):
        plt.show()

    def save(self, name="geodesic.png"):
        plt.savefig(name)
