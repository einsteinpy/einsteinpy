import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

from einsteinpy.metric import Schwarzschild


class StaticGeodesicPlotter:
    """
    Class for plotting static matplotlib plots.
    """

    def __init__(self, mass, time=0):
        self.mass = mass
        self.time = time

    def plot_attractor(self):
        plt.scatter(0, 0, color="black")

    def plot(self, pos_vec, vel_vec, end_lambda=10, step_size=1e-3):
        swc = Schwarzschild.from_spherical(pos_vec, vel_vec, self.time, self.mass)

        vals = swc.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]

        time = np.array([coord[0] for coord in vals])
        r = np.array([coord[1] for coord in vals])
        theta = np.array([coord[2] for coord in vals])
        phi = np.array([coord[3] for coord in vals])

        pos_x = r * np.cos(phi)
        pos_y = r * np.sin(phi)

        plt.scatter(pos_x, pos_y, s=1, c=time, cmap="Oranges")

    def show(self):
        plt.show()
