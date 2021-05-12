#requires transmit or reflect
def boundaryDetection(old_position, new_position,matrix):
    if(matrix[old_position.x,old_position.y,old_position.z].tissue==matrix[new_position.x,new_position.y,new_position.tissue.z]):
        return
    else:
        #transmit_or_reflect()//insert the function here when it's done
        return
def rouletteSurvive(photon, treshHold):
    """Determines if a photon survives or not based on its weight and a set reshold
    /aight Cosmin here: Please check the method I made in the Photon class. I think we should talk about it."""
    if photon.weight < tresHold:
        a = rnd.randint(0, 1/treshHold) #This one is very breakable. Must be improved
        if a == 1:
            photon.weight = 1
    else:
        return