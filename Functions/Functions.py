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

def euclid(x, y):
    while(y):
        x, y = y, x % y
  
    return x

def ReflectOrTransmit(n1, n2, theta1, Pincident):
    """
    Takes incident angle in radians, refractive inedecies of bodies and incident power. 
    Returns reflected and transmited power along with refraction angle
    """
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
def boundaryDetection(old_position, new_position,matrix,Pincident):
    #assuming tags are inside the matrix
    #Problem: what if the boundary is in between two identical points?
    '''
    alternative boundary detection: (use this only if the other one doesn't yield accurate results and you can lend processing time)
    sem=false
    directionVector.x=new_position.x-old_position.x
    directionVector.y=new_position.y-old_position.y
    directionVector.z=new_position.z-old_position.z
    k=euclid(directionVector.x,euclid(directionVector.y,directionVector.z))
    factor=0
    while factor*k<directionVector.x and factor*k<directionVector.y and factor*k<directionVector.z:
        for i in range(k):
            for j in range(k):
                for l in range(k):
                    if(matrix[old_position.x,old_position.y,old_position.z]==matrix[old_position.x+k*factor+i,old_position.y+k*factor+j,old_position.z+k*factor+z]):
                        sem=true
        factor+=1
    
    '''
    if(matrix[old_position.x,old_position.y,old_position.z]==matrix[new_position.x,new_position.y,new_position.tissue.z]):
        return new_position,0,Pincident
    else:
        rho=(new_position.x**2+new_position.y**2+new_position.z**2)**0.5
        #toDo: calculate incident angle, use matrix to get n1, n2, and indicent power? 
        #I'll just assume the incident angle is in radians.
        phi=math.atan((((new_position.x-old_position.x)**2+(new_position.y-old_position.y)**2)**0.5)/(new_position.z-old_position.z))
        theta=math.atan((new_position.y-old_position.y)/(new_position.x-old_position.x))
        #arrange call with Carlos and Sigurd, Get n1,n2, incident power from matrix.
        #boundary angles
        #also assumed phi is unchanged
        #verify that no angles are actuall cos/sin
        #double check coordinate systems
        theta2, PRefT, PTraT=ReflectOrTransmit(n1,n2,theta,Pincident)
        phi2, PRefP, PTraP=ReflectOrTransmit(n1,n2,phi,Pincident)
        #I am really not sure of how I am combining those two
        PRef=((PRefT+PRefP)/(2*Pincident))*Pincident
        PTra=((PTraT+PTraP)/(2*Pincident))*Pincident
        #converting from spherical to cartesian coordinates
        new_position.x=rho*math.sin(phi2)*math.cos(theta2)
        new_position.y=rho*math.sin(phi2)*math.sin(theta2)
        new_position.z=rho*math.cos(phi2)
        #returning data
        return new_position, PRef, PTra
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