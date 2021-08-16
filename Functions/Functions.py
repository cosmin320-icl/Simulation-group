#requires transmit or reflect
def boundaryDetection(old_position, new_position,matrix):
    if(matrix[old_position.x,old_position.y,old_position.z].tissue==matrix[new_position.x,new_position.y,new_position.tissue.z]):
        return
    else:
        #transmit_or_reflect()//insert the function here when it's done
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