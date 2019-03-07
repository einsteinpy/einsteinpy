import random

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from einsteinpy.metric import Schwarzschild
from einsteinpy.utils import schwarzschild_radius


class ScatterGeodesicPlotter:
    """
    Class for plotting static matplotlib plots.
    """

    def __init__(self, mass, time=0 * u.s):
        self.mass = mass
        self.time = time
        self._attractor_present = False

    def _plot_attractor(self):
        self._attractor_present = True
        plt.scatter(0, 0, color="black")

    def plot(self, pos_vec, vel_vec, end_lambda=10, step_size=1e-3):
        swc = Schwarzschild.from_spherical(pos_vec, vel_vec, self.time, self.mass)

        vals = swc.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]

        time = vals[:, 0]
        r = vals[:, 1]
        # Currently not being used (might be useful in future)
        # theta = vals[:, 2]
        phi = vals[:, 3]

        pos_x = r * np.cos(phi)
        pos_y = r * np.sin(phi)

        plt.scatter(pos_x, pos_y, s=1, c=time, cmap="Oranges")

        if not self._attractor_present:
            self._plot_attractor()

    def show(self):
        plt.show()

    def save(self, name="scatter_geodesic.png"):
        plt.savefig(name)


class StaticGeodesicPlotter:
    """
    Class for plotting static matplotlib plots
    """

    def __init__(self, mass, time=0 * u.s, ax=None):
        self.ax = ax
        if not self.ax:
            _, self.ax = plt.subplots(figsize=(6, 6))
        self.time = time
        self.mass = mass
        self._attractor_present = False

    def plot_trajectory(self, pos_vec, vel_vec, end_lambda, step_size, color):
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

        lines = self.ax.plot(x, y, "--", color=color)

        return lines, x[-1], y[-1]

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

        lines, x0, y0 = self.plot_trajectory(
            pos_vec, vel_vec, end_lambda, step_size, color
        )

        l, = self.ax.plot(x0, y0, "o", mew=0, color=lines[0].get_color())
        lines.append(l)

        self.ax.set_xlabel("$x$ (km)")
        self.ax.set_ylabel("$y$ (km)")
        self.ax.set_aspect(1)

        return lines

    def save(self, name="static_geodesic.png"):
        plt.savefig(name)
