import numpy as np

#daaph7 from BoostddG, seven different properties and difference is calculated
def get_daaph7(deleted,induced):
    amino_dict = {
        'A': [-0.35, -0.68, -0.677, -0.171, -0.17, 0.9, -0.476],
        'C': [-0.14, -0.329, -0.359, 0.508, -0.114, -0.652, 0.476],
        'D': [-0.213, -0.417, -0.281, -0.767, -0.9, -0.155, -0.635],
        'E': [-0.23, -0.241, -0.058, -0.696, -0.868, 0.9, -0.582],
        'F': [0.363, 0.373, 0.412, 0.646, -0.272, 0.155, 0.318],
        'G': [-0.9, -0.9, -0.9, -0.342, -0.179, -0.9, -0.9],
        'H': [0.384, 0.11, 0.138, -0.271, 0.195, -0.031, -0.106],
        'I': [0.9, -0.066, -0.009, 0.652, -0.186, 0.155, 0.688],
        'K': [-0.088, 0.066, 0.163, -0.889, 0.727, 0.279, -0.265],
        'L': [0.213, -0.066, -0.009, 0.596, -0.186, 0.714, -0.053],
        'M': [0.11, 0.066, 0.087, 0.337, -0.262, 0.652, -0.001],
        'N': [-0.213, -0.329, -0.243, -0.674, -0.075, -0.403, -0.529],
        'P': [0.247, -0.9, -0.294, 0.055, -0.01, -0.9, 0.106],
        'Q': [-0.23, -0.11, -0.02, -0.464, -0.276, 0.528, -0.371],
        'R': [0.105, 0.373, 0.466, -0.9, 0.9, 0.528, -0.371],
        'S': [-0.337, -0.637, -0.544, -0.364, -0.265, -0.466, -0.212],
        'T': [0.402, -0.417, -0.321, -0.199, -0.288, -0.403, 0.212],
        'V': [0.677, -0.285, -0.232, 0.331, -0.191, -0.031, 0.9],
        'W': [0.479, 0.9, 0.9, 0.9, -0.209, 0.279, 0.529],
        'Y': [0.363, 0.417, 0.541, 0.188, -0.274, -0.155, 0.476],
    }
    aapmt = np.array(list(amino_dict[induced]),dtype=np.float)
    aapwt = np.array(list(amino_dict[deleted]),dtype=np.float)
    daaph7 = aapmt-aapwt
    return list(daaph7)