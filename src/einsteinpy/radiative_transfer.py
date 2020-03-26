import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.coordinates import SphericalConversion
from einsteinpy.integrators import RK45

'''create a class photon which would include a position vector and a velocity vector'''
class radiative_transfer:
    '''This class holds methods to compute pitch angle , emmission coeffecient and radiative integration '''
    def __init__(self,photons,metric): 
        self.photon_pos=photon.position_vector #list of postions of individual photons
        self.photon_velo=photon.velocity_vector #list of velocities of individual photons
        self.photon_freq=photon.frequency #list of photon frequencies in plasma frame
        
    def pitch_angle(self):
        float p_ang = 0.0
        '''
        Determine pitch angle input params
        Variables holding sin(2*p_ang) and cos(2*p_ang) obtained using dot product between 
        magnetic field and polarisation basis
        Finally obatin p_ang from the previous results
        '''
        return p_ang 

    def emmission_coeffecient(self): #more parameters to be included
        float e_coeff = 0.0 #will hold the resultant emmission coeffecient
        '''
        e_c=Charge of electron
        c=Speed of light
        e_m=Electron Mass
        p_c=Planck Constant
        '''
        
        return e_coeff   
    
    def rad_integrate(self):
        '''
        This method involves defining the radiative transfer eqn integral and obtaining its solution
        dI/dlambda-define transfer eqn 
        Stokes parameters-I,Q,U,V
        Polarised Emissivities-e_j,e_q,e_u,e_v
        Absorption Coeffecients-a_i,a_q,a_u,a_v
        Define a method to obtain Emissivity and Intensity Vectors
        '''
    def rad_trans(self,dlambda,icur,intensity,tau):
        float dlambda_f=[dlambda/freq for dlambda,freq in zip(dlambda,self.photon_freq)] #frequency normalised values
        float e_c=emmission_coeffecient()








