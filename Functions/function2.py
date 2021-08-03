import random as rando
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
        delta_weight = (coeff_absorb / coeff_scatter) * self.weight
        #return self.weight - delta_weight (should the method go this far?)
        return delta_weight

    def scatter ():
        """"""
        pass