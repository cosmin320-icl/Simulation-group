import numpy as np
from numpy import random as rnd
import math

class PhotonClass():
    #from PhotonClass.py

    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
        return

    def removeweight (self, coeff_absorb, coeff_scatter):
        """inputs: coefficient of absorption, coefficient of scattering.
        Returns value (to be used as deltaweight),
        which is the amount of energy lost due to scattering at the end of a step."""
        delta_weight = (coeff_absorb / (coeff_absorb + coeff_scatter)) * self.weight

        #Update Luminosity 
        #Voxels[self.position[0] , self.position[1] , self.position[2] , 1] += delta_weight

        self.weight -= delta_weight
        return delta_weight


    def scatter (self, g):
        """
        g: scattering anisotropy (-1≦g≦1). If g is 0, this generally indicates that 
        the scattering is isotropic (i.e. no weighting in the scattering direction).
        If g approaches 1, this indicates that the scattering is primarily in the forward direction. (Wikipedia)
        ***
        updates direction"""

        # angles theta and phi
        if g == 0 :
            theta = math.acos(1 - (2 * rnd.rand()))
        else:
            theta = math.acos((1 / (2 * g)) * (1 + g**2 - ((1 - (g**2)) / (1 - g + 2 * g * rnd.rand()))))

        phi = 2 * math.pi * rnd.rand()

        #new directions (cannot update directions immediately as old values required to calculate new direction)
        new_myu_x = (math.sin(theta) * (self.direction[0] * self.direction[2] * math.cos(phi) - self.direction[1] * math.sin(phi))) / ((1 - (self.direction[2]**2))**0.5) + (self.direction[0] * math.cos(theta))
        new_myu_y = (math.sin(theta) * (self.direction[1] * self.direction[2] * math.cos(phi) - self.direction[0] * math.sin(phi))) / ((1 - (self.direction[2]**2))**0.5) + (self.direction[1] * math.cos(theta))
        new_myu_z = (-1) * (1 - self.direction[2]**2)**0.5 * math.sin(theta) * math.cos(phi) + self.direction[2] * math.cos(theta)
        
        #update direction
        self.direction[0] = new_myu_x
        self.direction[1] = new_myu_y
        self.direction[2] = new_myu_z

        # for special cases μ_z = 1 or -1? (as stated in wikipedia)
        return

def lum_update (Voxel_Matrix , photon , coeff_absorb, coeff_scatter):
    """updates the luminosity of matrix created by matrix maker with delta_weight (energy absorbed by cell)
    inputs: the matrix itself and the photon class
    
    **point to note** until further update, this function also does the weight subtraction i.e. absorption as well
        ↳in the future this function would probably be something like "absorption", which does the whole process for absorption (there'll then be another function for emission)"""

    dWeight = photon.removeweight(coeff_absorb, coeff_scatter)

    Voxel_Matrix[photon.position[0] , photon.position[1] , photon.position[2] ,1] += dWeight
    #4th index is always 1, as you're always manipulating the second 3D matrix which corresponds to luminosity
    
    return


Voxels = np.zeros((2,2,3,2), dtype = float) #assign_tags(sMatrix, materials) ←retrieves tag matrix from MatrixMaker.py in IO Group
#print(Voxels)
wan = PhotonClass (1, [1,1,1], [2,3,4])

lum_update(Voxels, wan, .5, .5)

print(Voxels)