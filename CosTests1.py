#Assume source location matrix looks like this
SLM = [{"x":1,"y":2,"z":3},{"x":4,"y":5,"z":6}] 

 

#def SourceDirectionCosine
DirectionCosineArray=[]
XYZ=["x","y","z"]
for elements in SLM:
    rmsOfCoordiates = (elements["x"]**2+elements["y"]**2+elements["z"]**2)**0.5
    DirectionCosineXYZ = []
    for coordinate in XYZ:
        DirectionCosineXYZ.append(elements[coordinate]/rmsOfCoordiates)
    DirectionCosineArray.append(DirectionCosineXYZ)
#Return DirectionCosineArray
a=input("bruh?")
print(a)