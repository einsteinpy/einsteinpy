import numpy as np
import scipy as si
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import integrate

# ds^2 = -(1-2*M/r)dt^2 + 1/(1-2*M/r)dr^2 + r^2*(d\theta^2 + sin(\theta)^2*d\phi^2)
M = 1

# Number of rays (n_rays) with the impact Paramter B
n_rays = 5000  # No. of rays
bcrit = 3 * np.sqrt(3) * M  # Critical Impact Paramter
bcrit = 5.196152422706632
obs = 30  # Observers Distance
B = np.linspace(bcrit, obs, n_rays)
# Solving for the turning points (r_tp) using Relation between the turning point and the impact parameter for the null geodesics [r_tp/(1-2*M/r) = b]. Only those light rays will be having turning points which are above the critical impact paramter (bcrit) which is 3*(3)**0.5*M, and all the rays below that falls into the spacetime. Its is polynomial equation so the turning point will the real maximum root (having maximum effective potential).

bn = np.zeros(len(B))
R = []
b_fin = []
for i in B:
    M = 1
    eqn1 = lambda r: r / ((1 - (2 * M / r))) ** 0.5 - i
    root = si.optimize.newton(eqn1, 0.1)
    if np.isreal(root) == True:
        R.append(np.real(root))
        b_fin.append(i)

z = np.asarray(list(zip(b_fin, R)))


# Defining the function for emissivity (Intensity).

# Using the backward raytracing the rays will shot back in time from the observer till the source and then the intensity for the will be calculated on the basis of red or blue shift which will chane a positive/negative sign in the equation.
def intensity_blue_sch(r, b):
    GTT = (1 - (2 * M / r))
    GRR = (1 - (2 * M / r)) ** (-1)

    KRKText = ((GTT / GRR) * (1 - (b ** 2 * GTT / (r ** 2)))) ** 0.5
    Gblue = ((1 / GTT) - KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5) ** (-1)
    Iblue = -(Gblue ** 3) * (GTT / (r ** 2)) * (1 / KRKText)
    return Iblue


def intensity_red_sch(r, b):
    GTT = (1 - (2 * M / r))
    GRR = (1 - (2 * M / r)) ** (-1)

    KRKText = ((GTT / GRR) * (1 - (b ** 2 * GTT / (r ** 2)))) ** 0.5
    Gred = ((1 / GTT) + KRKText * (GRR / GTT) * ((1 - GTT) / (GTT * GRR)) ** 0.5) ** (-1)
    Ired = (Gred ** 3) * (GTT / (r ** 2)) * (1 / KRKText)
    return Ired


#source and observer are considered to be on same distances from the centre
#Numerically integrating the intensity along the null geodesics from the observer till the source.

Int1 = []
for i in np.arange(len(z)):
    b = z[i,0]
    val1, err1 = si.integrate.quadrature(intensity_blue_sch,obs,z[i,1],args=(b,))
    val2, err2 = si.integrate.quadrature(intensity_red_sch,z[i,1],obs,args=(b,))
    val_fin1 = val1+val2
    Int1.append(val_fin1)

# Integration the inensity from the singularity till the horizon (It will only contribute to the redshift part.)

sing = 0.01  # Singularity
B2 = np.linspace(sing, bcrit, len(b_fin))
k1 = []

# Intensities from horzion to the observer.
for i in B2:
    b = i
    val3, err3 = si.integrate.quadrature(intensity_red_sch, 2, 30, args=(b,))
    k1.append(val3)

Int2 = k1
Int3 = Int2+list(Int1)

fb1 = list(B2) + list(b_fin)
fb2 = list(-np.asarray(B2)) + list(-np.asarray(b_fin)) #Just to make the plot symmetric on -x axis
plt.plot(fb1,Int3,'r')
plt.plot(fb2,Int3,'r')
plt.xlabel('Impact Paramter (b)')
plt.ylabel('Intensity (Emissivity)')
plt.title('Intensity Plot')
plt.show()

#Making the Shadow plot (Polar)
# radius and theta
theta1 = np.linspace(0,2*np.pi,len(fb1))
r1, theta1 = np.meshgrid(fb1, theta1)

#Intensities
values1, values2 = np.meshgrid(Int3,Int3)

xx= r1*np.cos(theta1)
yy = r1*np.sin(theta1)

plt.pcolormesh(xx,yy,values1,cmap=plt.cm.hot,shading='gouraud')
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.gca().set_aspect('equal', adjustable='box')
