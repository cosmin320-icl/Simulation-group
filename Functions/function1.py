import math
import numpy as np
from numpy import random as rnd

"""
Function 1([inputs]): 
#s=Determine Step Size 
    s=-lnE/ut

PhotonClass.MovePhoton(s): 
Update photon position 
Return PhotonClass """

class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
        return

    def MovePhoton(self, step):
        #point to consider: what if direction vector wasn't a unit vector?
        #or are we limiting to unit vectors?
        self.newPosition=np.array(self.position)+(np.array(self.direction)*step)
        return self.newPosition

    #Function 1 but as a class method
    def m_Function1 (self, cAbsorption, cScattering):
        cInteraction=cAbsorption+cScattering
        s=(-(math.log(rnd.random())/cInteraction))
        print(s)
        self.position=self.MovePhoton(s)
        return self


def Function1Test(cInteraction):
    testPhoton=PhotonClass(1,[0,0,0],[0,0,1])
    testPhoton=testPhoton.m_Function1(cInteraction)
    assert testPhoton.position[0]==testPhoton.position[1]==0
    return

    

#how do you update a class attribute from outside the class?
#also recall documentation for random library
# def Function1 (cInteraction, photon):
#     s=(-(math.log(rnd.random())/cInteraction))
#     photon.position=photon.MovePhoton(s)
#     return photon

# test=PhotonClass(1,[0,0,0],[0,0,1])
# test.m_Function1(10)
# print(test.position)

Function1Test(10)

#Function1Test(10)

class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
        return
    
#Assume source location matrix looks like this
SLM = [{"x":1,"y":2,"z":3},{"x":4,"y":5,"z":6}] 

#Tests for Taiga /Cosmin
def Function1Test(Function1):
    testPhoton=PhotonClass(1,[0,0,0],[0,0,1])
    testPhoton=Function1(PhotonClass)
    assert testPhoton.position[0]==testPhoton.position[1]==0
    return
