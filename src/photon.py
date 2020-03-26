import astropy.units as u
import numpy as np

from einsteinpy import constant

''' any ideas for any existing base class this can inherit ?'''

class photon:

    def __init__(self,vec): #initialisation of class variables from a single vector
        self.frequency=np.zeros(shape=vec.shape,dtype=vec.dtype)
        self.position_vector=np.zeros(shape=vec.shape,dtype=vec.dtype)
        self.velocity_vector=np.zeros(shape=vec.shape,dtype=vec.dtype)
    def input_photon_params(self):
        '''define input type and required dimensions'''
