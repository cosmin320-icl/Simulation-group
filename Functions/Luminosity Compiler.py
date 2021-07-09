import numpy as np
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