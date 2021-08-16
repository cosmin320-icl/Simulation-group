import random as rand
import numpy as np
import math

class PhotonClass():
    def __init__(self,weight,position,direction):
        """photon direction must be a 3D unit vector, therefore the magnitude has to be 1 
        (any greater, and the program starts to spit complex numbers)
        """
        self.weight=weight
        self.position=np.array(position)
        self.direction=np.array(direction)
        return

    def rouletteSurvive(self,treshold):
        if(self.weight<treshold):
            #m is given by the 1 in m chances to survive. What value should m have?
            m=5
            if(rando.random()<(1/m)):
                self.weight*=m
            else:
                self.weight=0
            return
        else:
            return

    def MovePhoton(self, coeff_absorb, coeff_scatter):
        #point to consider: what if direction vector wasn't a unit vector?
        #or are we limiting to unit vectors?
        step=(-(math.log(rand.random())/(coeff_absorb + coeff_scatter)))
        self.position = (self.position) + (self.direction)*step
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

