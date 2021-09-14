import numpy as np
from numpy import random as rnd

def pointOnPlane(vector, point, xLength = 1, yLength = 1, circular = False):
    "inputs 1*3 arrays for point and vector, to generate a point within the plane these define and the constraints xLength and yLength"
    defaultPlane = np.array([0, 0, 1]) #plane where coordinates are generated
    vector = vector / sum(vector**2) ** (1/2) #normalises input vector

    #generate an x and y coordinate on the (0,0,1) plane going thorugh the origin
    if circular != False:
        #to generate points on plane within a circle
        x = rnd.random()*circular - circular/2
        circular2 = 2*((circular/2)**2 - x**2)**(1/2)
        y = rnd.random()*circular2 - circular2/2
    else:
        #to generate points on plane within a rectangle
        x = rnd.random()*xLength - xLength/2
        y = rnd.random()*yLength - yLength/2
    
    #finds the angle between the plane 
    angle = -np.arccos(np.dot(defaultPlane, vector))
    #small quirk that angle is negative, I think it is because of the order of the cross-product
    
    if angle == 0:
        return np.array([x + point[0], y + point[1], 0 + point[2]]) #clause to avoid errors
        
    #finds line of rotation
    r = np.cross(vector, defaultPlane)
    #normalises this vector
    r = r / sum(r ** 2) ** (1/2)
    

    

    #rotates to wanted points and adjusts for coordinates
    xr = point[0] + ((1-np.cos(angle))*r[0]**2 + np.cos(angle))*x + ((1-np.cos(angle))*r[0]*r[1] - np.sin(angle)*r[2])*y 
    yr = point[1] + ((1-np.cos(angle))*r[0]*r[1] + np.sin(angle)*r[2])*x + ((1-np.cos(angle))*r[1]**2 + np.cos(angle))*y  
    zr = point[2] + ((1-np.cos(angle))*r[0]*r[2] - np.sin(angle)*r[1])*x + ((1-np.cos(angle))*r[1]*r[2] + np.sin(angle)*r[0])*y
    
    #Organises coordinates xr, yr, zr into a point on the wanted plane
    onPlane = np.array([xr, yr, zr])
    
    #not done, will add to integrate with photonclass
    return onPlane

def determineStepSize(u_t, scaleFactor):
    "u_t is the sum of scattering and absorption coefficient of a voxel"
    #Assuming the random number is between 0 and 1
    stepSize = - np.log(rnd.random()) / u_t 
    #Stepsize is in centimeters, scaled to fit with matrix parameters
    stepSize *= scaleFactor
    return stepSize

import random as rando
class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        #the direction vector gets normalised
        self.direction = direction / (sum(direction**2))**(1/2)
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
        totalInteractionCoeff = absoprtionCoeff + scatteringCoeff
        deltaW = (absoprtionCoeff / totalInteractionCoeff) * self.weight
        self.weight-=deltaW
        return deltaW

    def movePhoton(self, stepSize):
        self.position = self.position + self.direction * stepSize
        return 

    def scatterPhoton(self, g = 0.9):
        """g is the anisotropy, assumed to be 0.9 for all biological tissues based on data we found"""
        #Normalising self.direction to be on the safe side
        self.direction = self.direction / (sum(self.direction**2))**(1/2)
        cosTheta = 1.0/(2.0*g) * ((1.0 + g**2) - ((1.0 - g**2.0)/(1 - g + 2*g*rnd.random()))**2)
        sinTheta = np.sin(np.arccos(cosTheta))
        pAngle = 2*np.pi*rnd.random()

        if self.direction[2]== 1:
            self.direction[0] = sinTheta * np.cos(pAngle)
            self.direction[1] = sinTheta * np.sin(pAngle)
            self.direction[2] = cosTheta

        elif self.direction[2] == -1:
            self.direction[0] = sinTheta * np.cos(pAngle)
            self.direction[1] = -sinTheta * np.sin(pAngle)
            self.direction[2] = -cosTheta    
        else:
            self.direction[0] = sinTheta *(self.direction[0]*self.direction[2]*np.cos(pAngle) - self.direction[1]*np.sin(pAngle)) / (1 -  self.direction[2]**2)**(1/2) + self.direction[0]*cosTheta
            self.direction[1] = sinTheta *(self.direction[1]*self.direction[2]*np.cos(pAngle) - self.direction[0]*np.sin(pAngle)) / (1 -  self.direction[2]**2)**(1/2) + self.direction[1]*cosTheta   
            self.direction[2] = - (1 - self.direction[2]**2)**(1/2)*sinTheta*np.cos(pAngle) + self.direction[2] * cosTheta
            
        #Normalising self.direction to be on the safe side
        self.direction = self.direction / (sum(self.direction**2))**(1/2)
        return self.direction
