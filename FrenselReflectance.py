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