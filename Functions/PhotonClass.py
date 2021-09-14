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
    def removeWeight(self, absoprtionCoeff, scatteringCoeff):
        totalInteractionCoeff=absoprtionCoeff+scatteringCoeff
        deltaW=(absoprtionCoeff/totalInteractionCoeff)*self.weight
        self.weight-=deltaW
        return deltaW
pho=PhotonClass(100,[1,2,1],[1,1,4])
deltaW=pho.removeWeight(0.2,0.22)
print(pho.weight+deltaW)