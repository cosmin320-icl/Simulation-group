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

        #return self.weight - delta_weight (should the method go this far?)
        return delta_weight

        """Matrix[position].luminosity+=deltaWeight 

        Photon.Class.scatter() 

        Return LuminosityQueue.put(Update Luminosity in List),[rework: position,luminosity], photonClass """


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

wan = PhotonClass (1, [1,1,1], [2,3,4])
wan.scatter(0.5)
print(wan.direction)
#output : [(1.055866209337768-1.4880055503868064j), (1.583799314006652-1.7672229971970363j), (2.111732418675536-1.7442434181131221j)]
#???

# update luminosity?? what is the matrix?