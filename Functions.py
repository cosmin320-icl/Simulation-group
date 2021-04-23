#requires transmit or reflect
def boundary detection(old_position, new_position,matrix):
    if(matrix[old_position.x,old_position.y,old_position.z].tissue==matrix[new_position.x,new_position.y,new_position.tissue.z]):
        return
    else:
        transmit_or_reflect()
        return
