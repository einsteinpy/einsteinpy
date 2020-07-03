import numpy as np
from matplotlib import pyplot as plt


class ShadowPlotter:
    """
    Class for plotting and visualising shadows
    """

    def __init__(self, shadow, is_line_plot=True):
        """
        Constructor for plotter.

        Parameters
        ----------
        shadow : ~einsteinpy.rays.Shadow
            The shadow object
        is_line_plot : bool, optional
            If the plot is a line plot or a contour plot. Defaults to True.
        """
        self.shadow = shadow
        self.is_intensity_plot = is_line_plot

    def plot(self):
        """
        Plots the shadow.
        """
        if self.is_intensity_plot:
            plt.plot(self.shadow.fb1, self.shadow.intensity, "r")
            plt.plot(self.shadow.fb2, self.shadow.intensity, "r")
            plt.xlabel("Impact Paramter (b)")
            plt.ylabel("Intensity (Emissivity)")
            plt.title("Intensity Plot")
        else:
            theta1 = np.linspace(0, 2 * np.pi, len(self.shadow.fb1))
            self.r1, self.theta1 = np.meshgrid(self.shadow.fb1, theta1)
            self.values1, self.values2 = np.meshgrid(
                self.shadow.intensity, self.shadow.intensity
            )

    def show(self):
        """
        Shows the plot.
        """
        if self.is_intensity_plot:
            plt.show()
        else:
            xx = self.r1 * np.cos(self.theta1)
            yy = self.r1 * np.sin(self.theta1)
            plt.figure(figsize=(7, 7))
            plt.pcolormesh(xx, yy, self.values1, cmap=plt.cm.afmhot, shading="gouraud")
            plt.title("Schwarzschild Black Hole")
            plt.gca().set_aspect("equal", adjustable="box")
