from numpy import random as rnd
import numpy as np

def function0(vector, point, xLength = 1, yLength = 1, circular = False):
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