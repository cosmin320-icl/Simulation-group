#Ahmed's Stuff
import numpy as np
import numpy.random as random

# set dimensions of test matrix
x= 5
y= 5
z= 3

photonArray = [] # assume photon objects are stores in a list
photonsPerNano = 10

nanoparticleFrequency = 10 # how many particles there are

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

# initiate arrays
Matrix = CreateCoordinateArray(x, y, z)
nanoparticleArray = CreateParticleArray(nanoparticleFrequency, Matrix)

def updateValues(Matrix, nanoparticleArray):
    for i in nanoparticleArray:
    # access voxel, then do...?
    Matrix[i[0]][i[1]][i[2]] = 1 # testing if access of matrix at positions are correct
    # create photon using PhotonClass(). photon objects stored in list?
    # for i in nanoparticleArray: // i is the position
    # for j in range(photonsPerNano)
    # photonArray.append(PhotonClass(i, weight??, direction?? ))) // weight and direction are beyond my understanding
    return Matrix


print("Matrix before update")
print(Matrix)

Matrix = updateValues(Matrix, nanoparticleArray)

print("Matrix after update")
print(Matrix)