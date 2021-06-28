#requires transmit or reflect
#by Masaki

import math
### Material information (Dummy values) 
n1=1.000293     #Refractive index (First medium)
n2=1.333     #Refractive index (Second medium)
theta1 = 0.524     #radian
Pincident = 1     #Power(or energy?) of incident light
 
#Fresnel_equations
Rs = ((n1*math.cos(theta1)-n2*math.sqrt(1-(n1/n2*math.sin(theta1))**2))/(n1*math.cos(theta1)+n2*math.sqrt(1-(n1/n2*math.sin(theta1))**2)))**2
Rp = ((-n2*math.cos(theta1)+n1*math.sqrt(1-(n1/n2*math.sin(theta1))**2))/(n2*math.cos(theta1)+n1*math.sqrt(1-(n1/n2*math.sin(theta1))**2)))**2

 

Reff = 0.5*(Rs+Rp)     #Effective reflectance
T = 1-Reff     #Transmittance
PRef = Pincident*Reff     #Reflected Power = incident P * R
PTra = Pincident*T     #Transmitted P = incident P * T
print("Effective reflectance: {}, Transmittance: {}".format(Reff,T))
print("Power after reflection: {}, P after Transmission: {}".format(PRef,PTra))

 

#New direction
theta2 = math.asin(n1/n2*math.sin(theta1))
print("New direction, theta2:", theta2)

import math 

def ReflectOrTransmit(n1, n2, theta1, Pincident):
    """Takes incident angle in radians, refractive inedecies of bodies and incident power. 
    Returns reflected and transmited power along with refraction angle"""
    #tried to make a function of the stuff above
    Rs = ((n1*math.cos(theta1)-n2*math.sqrt(1-(n1/n2*math.sin(theta1))**2))/(n1*math.cos(theta1)+n2*math.sqrt(1-(n1/n2*math.sin(theta1))**2)))**2
    Rp = ((-n2*math.cos(theta1)+n1*math.sqrt(1-(n1/n2*math.sin(theta1))**2))/(n2*math.cos(theta1)+n1*math.sqrt(1-(n1/n2*math.sin(theta1))**2)))**2

    Reff = 0.5*(Rs+Rp)     #Effective reflectance
    T = 1-Reff     #Transmittance
    PRef = Pincident*Reff     #Reflected Power = incident P * R
    PTra = Pincident*T     #Transmitted P = incident P * T
    
    #new direction
    theta2 = math.asin(n1/n2*math.sin(theta1))
    
    return theta2, PRef, PTra
def boundaryDetection(old_position, new_position,matrix):
    #assuming tags are inside the matrix
    if(matrix[old_position.x,old_position.y,old_position.z]==matrix[new_position.x,new_position.y,new_position.tissue.z]):
        return
    else:
        #toDo: calculate incident angle, use matrix to get n1, n2, and indicent power? 
        #I'll just assume the incident angle is in radians.
        phi=math.atan((((new_position.x-old_position.x)**2+(new_position.y-old_position.y)**2)**0.5)/(new_position.z-old_position.z))
        theta=math.atan((new_position.y-old_position.y)/(new_position.x-old_position.x))
        #arrange call with Carlos and Sigurd, Get n1,n2, incident power from matrix.
        theta2, PRef, PTra=ReflectOrTransmit(n1,n2,phi,Pincident);
        #convert from polar to cartesian 
        #add old_coordinates
        
        return
        '''
def rouletteSurvive(photon, treshHold):
    """Determines if a photon survives or not based on its weight and a set reshold
    /aight Cosmin here: Please check the method I made in the Photon class. I think we should talk about it."""
    if photon.weight < tresHold:
        a = rnd.randint(0, 1/treshHold) #This one is very breakable. Must be improved
        if a == 1:
            photon.weight = 1
    else:
        return
        '''