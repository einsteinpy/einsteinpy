import numpy as np
import astropy.units as u
from einsteinpy.metric import Schwarzschild
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from einsteinpy import constant


M = 1.989e30 * u.kg
pos_vec = np.array([150e9, np.pi/2, 0.])
real_omega = 29951.68/150e9
vel_vec = np.array([0, 0., real_omega])


time = 0 * u.s
cl = Schwarzschild.from_values(pos_vec, vel_vec, time, M)
ANS = cl.calculate_trajectory(end_lambda=31558464*constant.c.value, OdeMethodKwargs={'stepsize':constant.c.value*3000})
# ANS = cl.calculate_trajectory(end_lambda=300000, OdeMethodKwargs={'vectorized':True})
ans = ANS[1]
print('calculate done')

r = np.array([t[1] for t in ans])
print(r[-1])

phi = np.array([t[3] for t in ans])
time = np.array([t[0] for t in ans])
x = r * np.cos(phi)
y = r*np.sin(phi)
plt.scatter(x,y, s=1, c=time, cmap='Oranges')
plt.scatter(0, 0 , color='black')
plt.show()