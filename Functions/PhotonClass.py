import random as rando
class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
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
    