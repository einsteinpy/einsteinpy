import numpy as np
import astropy.units as u
from einsteinpy.metric import Schwarzschild
# from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt


pos_vec = np.array([305.0, np.pi/2, 0.])
# vel_vec = np.array([0, 0, 6.6591005e-6])
vel_vec = np.array([0, 0, 6.6447e-6])
# vel_vec = np.array([0, 0, 0])
M = 5.972e25 * u.kg
time = 0 * u.s
cl = Schwarzschild.from_values(pos_vec, vel_vec, time, M)
ANS = cl.calculate_trajectory(end_lambda=300000, steplen =4.0)
ans = ANS[1]
print('calculate done')

r = np.array([t[1] for t in ans])
print(r[-1])

phi = np.array([t[3] for t in ans])
time = np.array([t[0] for t in ans])
x = r * np.cos(phi)
y = r*np.sin(phi)
plt.scatter(x,y, s=3, c=time, cmap='Oranges')
plt.scatter(0,0, color='black')
plt.show()
