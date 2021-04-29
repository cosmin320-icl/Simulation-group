class PhotonClass():
    def __init__(self,weight,position,direction):
        self.weight=weight
        self.position=position
        self.direction=direction
        return
    
#Assume source location matrix looks like this
SLM = [{"x":1,"y":2,"z":3},{"x":4,"y":5,"z":6}] 

def SourceDirectionCosine(SLM):
    DirectionCosineArray=[]
    XYZ=["x","y","z"]
    for elements in SLM:
        rmsOfCoordiates = (elements["x"]**2+elements["y"]**2+elements["z"]**2)**0.5
        DirectionCosineXYZ = []
        for coordinate in XYZ:
            DirectionCosineXYZ.append(elements[coordinate]/rmsOfCoordiates)
        DirectionCosineArray.append(DirectionCosineXYZ)
    return DirectionCosineArray
#I don't think I understand the code well enough. Ask Masaki for clarifications /Cosmin

#Tests for Taiga /Cosmin
def Function1Test(Function1):
    testPhoton=PhotonClass(1,[0,0,0],[0,0,1])
    testPhoton=Function1(PhotonClass)
    assert testPhoton.position[0]==testPhoton.position[1]==0
    return
#Tests for Taiga /Cosmin
def Function2ScatterTest(Function2Scatter,finalPosition,matrix, listOfLuminosities):
    testPhoton=PhotonClass(1,[0,0,0],[0,0,1])
    position, luminosity, testPhoton=Function2Scatter()
    desiredLuminosity=0
    desiredPositon=[0,0,0]
    semaL=True
    semaP=True
    #semaphore for luminosity:
    if(luminosity!=desiredLuminosity):
        semaL=False
    if(position!=[0,0,0]):
        semaP=False
    assert (semaL and semaP)==True
    return 
    #NOT YET DONE, REWORKS PROPOSED
