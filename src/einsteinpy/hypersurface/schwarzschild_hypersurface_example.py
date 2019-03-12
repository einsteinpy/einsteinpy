import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt

def gradient(mass, r):
    R = r / math.sqrt(1 - (2 * mass / r))
    num_one = (1 - (3 * mass / r))
    num_two = math.sqrt(((4 * mass * r - 9 * mass * mass) * R) / (r - 3 * mass) ** 2)
    deno = (math.sqrt(1 - (2 * mass / r)) ** 3)
    
    return num_one * num_two / deno


def R(mass, r):
    return r / math.sqrt(1 - (2 * mass / r))

def get_values(mass, alpha):
    x_axis = []
    y_axis = []
    r_initial = 3 * mass + 0.01 # just equal to schwarzschild radius
    r_step = mass / alpha
    
    # for values greater than schwarzschild radius
    z = 0
    r = r_initial
    while (r < 9 * mass):
        x_axis.append(R(mass, r))
        y_axis.append(z)
        z = z + gradient(mass, r) * r_step
        r = r + r_step
        
    # for values less than schwarzschild radius but greater than 9m/4
    z = 0
    r = r_initial
    while (r > (9 * mass / 4)):
        x_axis.append(R(mass, r))
        y_axis.append(z)
        z = z + gradient(mass, r) * r_step
        r = r - r_step
    
    return x_axis, y_axis


def get_values_surface(mass, alpha):
    r_initial = 3 * mass + 0.01 # just equal to schwarzschild radius
    r_step = mass / alpha
    phi_values = np.linspace(0, 2*np.pi, 60)
    R_values = []
    z_values = []
    
    # for values greater than schwarzschild radius
    z = 0
    r = r_initial
    while (r < 10 * mass):
        R_values.append(R(mass, r))
        z_values.append(z)
        z = z + gradient(mass, r) * r_step
        r = r + r_step
    
    """
    # for values less than schwarzschild radius but greater than 9m/4
    z = 0
    r = r_initial
    while (r > (9 * mass / 4)):
        R_values.append(R(mass, r))
        z_values.append(z)
        z = z + gradient(mass, r) * r_step
        r = r - r_step
    """
    R_values = np.array(R_values)
    R_values, phi_values = np.meshgrid(R_values, phi_values)
    
    X = R_values * np.cos(phi_values)
    Y = R_values * np.sin(phi_values)
    x_len = X.shape[0]
    y_len = X.shape[1]
    
    
    Z = np.array(z_values)
    z_values = np.array(z_values)
    for i in range(0, x_len - 1):
        Z = np.concatenate((Z, z_values), axis = 0)
    
    Z.reshape((x_len, y_len))
    return X, Y, Z



x_values, y_values = get_values(10, 100)
plt.scatter(x_values, y_values)
plt.show()

# plot the 3d hypersurface
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')
X, Y, Z = get_values_surface(10, 100)
shape_tuple = X.shape
Z = Z.reshape((shape_tuple[0], shape_tuple[1]))
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
 #               cmap='cubehelix', edgecolor='none')
ax.plot_wireframe(X, Y, Z, color='black')
