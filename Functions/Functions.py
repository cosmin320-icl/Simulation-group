#Declare required modules
from numpy.matrixlib.defmatrix import matrix
from Functions.function0 import function0
import random as rando
from numpy import random as rnd
import numpy as np
import math
import numpy.random as random
#Declare PhotonClass. This is the class which stores all the properties of a photon and operations related to it
class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
        #End of init
        return
    def rouletteSurvive(self,treshold):
        #Roulette survive is created in order to maintain the mathematical validity of the system. It is simply a function which terminates photons such that the stochastic analysis remains valid.
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
        #this is the functiont that removes the weight scattered by the photon and returns the difference between the current and previous weight
        totalInteractionCoeff=absoprtionCoeff+scatteringCoeff
        deltaW=(absoprtionCoeff/totalInteractionCoeff)*self.weight
        self.weight-=deltaW
        return deltaW
    
    def m_Move_Photon (self, cAbsorption, cScattering):
        #Function 1 but as a class method
        cInteraction=cAbsorption+cScattering
        s=(-(math.log(rnd.random())/cInteraction))
        print(s)
        self.position=self.MovePhoton(s)
        return self
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
'''
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
'''
def euclid(x, y):
    #Euclid's algorithm. please ignore
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
#could this be vectorized in some way


def LuminosityCompiler(Matrix, LuminosityList):
    "puts luminosity of voxels into layers txt file, there will be one txt file per layer"
    for a in range(length):
        layer = np.zeroes(height, width)
        for i in range(0, height):
            for j in range(0, width):
                layer[i][j] = voxel.luminosity
    
        layertxt = open("Layer_{0:05d}.txt".format(i), "w+")
        np.savetxt("Layer_{0:05d}.txt".format(i), layer, fmt="%s")
        Snapshot.close()
    
    return
def Generate_Photon(vector, point, xLength = 1, yLength = 1, circular = False):
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
def Move_Photon(Photon):
    #DO NOT USE. REPLACED BY PhotonClass.m_Function1()
    #Hopefully this is in Taiga's code. I couldn't find it
    #s=DetermineStepSize()
    #PhotonClass.movePhoton(s)
    Photon.m_Function1()
    return Photon

def Determine_Weight(Matrix,photon,absoprtionCoeff,scatteringCoeff,Voxel_Matrix, LuminosityList, old_position, new_position,Pincident):
    deltaWeight=photon.removeWeight(absoprtionCoeff,scatteringCoeff)
    #difference between matrix and voxel matrix??
    #update luminosity in matrix
    lum_update (Voxel_Matrix , photon , absoprtionCoeff, scatteringCoeff)
    #PhotonClass.scatter()
    photon.scatter()
    #update luminosity queue
    LuminosityCompiler(Matrix, LuminosityList)
    #additional checks: If boundary, if transmit or reflect
    new_position, PRef,Ptra=boundaryDetection(old_position,new_position,Matrix, Pincident)
    #what exactly do we do with PRef?
    Pincident=Ptra
    return new_position, Pincident, LuminosityList,Matrix,Voxel_Matrix

def Check_Weight(photon,treshold,tresholdSurvive):
    if photon.weight<treshold:
        photon.rouletteSurvive(tresholdSurvive)
    return
def checkOutside(position,matrixBoundaries):
    #check if position is outside, return true
    return False
#Ahmed's stuff
def updateValues(Matrix, nanoparticleArray):
    for i in nanoparticleArray:
    # access voxel, then do...?
    Matrix[i[0]][i[1]][i[2]] = 1 # testing if access of matrix at positions are correct
    # create photon using PhotonClass(). photon objects stored in list?
    # for i in nanoparticleArray: // i is the position
    # for j in range(photonsPerNano)
    # photonArray.append(PhotonClass(i, weight??, direction?? ))) // weight and direction are beyond my understanding
    return Matrix
def CreateCoordinateArray(x, y, z):
    # returns test matrix
    dimensions = (y, x, z)
    array1 = np.zeros(dimensions, dtype = int)
    return array1

def CreateParticleArray(frequency, Matrix):
    # returns array of random voxel position in which nanoparticle == true. Used for testing
    array = []
    xArray = random.randint(0, np.shape(Matrix)[1], size = frequency)
    yArray = random.randint(0, np.shape(Matrix)[0], size = frequency)
    zArray = random.randint(0, np.shape(Matrix)[2], size = frequency)
    for i in range(frequency):
    randomCoordinate = (yArray[i], xArray[i], zArray[i])
    array.append(randomCoordinate)
    return array