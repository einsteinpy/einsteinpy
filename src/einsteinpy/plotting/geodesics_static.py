import random

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from einsteinpy.metric import Schwarzschild
from einsteinpy.utils import schwarzschild_radius


class StaticGeodesicPlotter:
    """
    Class for plotting static matplotlib plots
    """

    def __init__(self, mass, time=0 * u.s, ax=None):
        self.ax = ax
        if not self.ax:
            self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.time = time
        self.mass = mass
        self._attractor_present = False

    def get_trajectory(self, pos_vec, vel_vec, end_lambda, step_size):
        swc = Schwarzschild.from_spherical(pos_vec, vel_vec, self.time, self.mass)

        vals = swc.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]

        time = np.array([coord[0] for coord in vals])
        r = np.array([coord[1] for coord in vals])
        theta = np.array([coord[2] for coord in vals])
        phi = np.array([coord[3] for coord in vals])

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        return x, y

    def plot_attractor(self):
        if not self._attractor_present:
            self._draw_attractor()

    def _draw_attractor(self, min_radius=10 * u.km):
        radius = max(schwarzschild_radius(self.mass) * 10000, min_radius.to(u.km))
        color = "#ffcc00"
        self.ax.add_patch(mpl.patches.Circle((0, 0), radius.value, lw=0, color=color))

    def plot(
        self,
        pos_vec,
        vel_vec,
        end_lambda=10,
        step_size=1e-3,
        color="#{:06x}".format(random.randint(0, 0xFFFFFF)),
    ):

        self.plot_attractor()
        self._attractor_present = True

        pos_x, pos_y = self.get_trajectory(pos_vec, vel_vec, end_lambda, step_size)

        lines = self.ax.plot(pos_x, pos_y, "--", color=color)
        l, = self.ax.plot(pos_x[-1], pos_y[-1], "o", mew=0, color=lines[0].get_color())
        lines.append(l)

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines

    def animate(
        self,
        pos_vec,
        vel_vec,
        end_lambda=10,
        step_size=1e-3,
        color="#{:06x}".format(random.randint(0, 0xFFFFFF)),
        interval=50,
    ):

        pos_x, pos_y = self.get_trajectory(pos_vec, vel_vec, end_lambda, step_size)
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

        ani = FuncAnimation(
            self.fig, _update, frames=frames, interval=interval, blit=True
        )
        plt.show()

    def save(self, name="static_geodesic.png"):
        plt.savefig(name)
