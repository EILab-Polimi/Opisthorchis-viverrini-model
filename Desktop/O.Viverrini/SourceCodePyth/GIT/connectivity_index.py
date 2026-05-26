import numpy as np
from build_connectivity_matrix import build_connectivity_matrix
def connectivity_index(OCNRes, nNodes, dam_pass, alpha):
    nres = OCNRes['nres']
    river = []
    for i in range(len(OCNRes['FD']['A'])):
        if OCNRes['FD']['A'][i] >= OCNRes['thrA']:
            down = int(OCNRes['FD']['downNode'][i]-1)  # 0-based
            dist = np.sqrt((OCNRes['FD']['X'][i] - OCNRes['FD']['X'][down])**2 +
                           (OCNRes['FD']['Y'][i] - OCNRes['FD']['Y'][down])**2) \
                   if i != OCNRes['outlet'] - 1 else 1
            river.append([i, down, dist, OCNRes['ReservoirsPixel'][i]])
    river = np.array(river)
    pixel2rivnode_full = OCNRes['MSC'].flatten()
    pixel2rivnode_full = pixel2rivnode_full.T
    LambdaF_mat, Lambda1 = build_connectivity_matrix(river, nNodes, dam_pass, alpha, pixel2rivnode_full,nres)
    return LambdaF_mat, Lambda1

