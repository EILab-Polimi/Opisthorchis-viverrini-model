import numpy as np
import scipy.io as sio
import os

def build_OCN_Reservoirs(name, thrA, nReservoirs, seed_value, habitatsuitability=None, seedchange = 0):
    """
    Traduzione Python 1:1 della funzione MATLAB build_OCN_Reservoirs.
    Ricrea la stessa struttura OCN con le stesse variabili e comportamento.
    """

    # Setting seed for reproducibility
    np.random.seed(seed_value)

    # load file .mat
    data = sio.loadmat(os.path.join("dataOCN", name))

    # Extracting variables from .mat
    X = data['X'].squeeze()
    Y = data['Y'].squeeze()
    CTC = data['CTC'].squeeze()
    SCX = data['SCX'].squeeze()
    SCY = data['SCY'].squeeze()
    outlet = data['outlet'].squeeze()
    FD_A = data['FD_A']
    FD_X = data['FD_X']
    FD_Y = data['FD_Y']
    FD_downNode = data['FD_downNode'].squeeze()
    downNode = data['downNode'].squeeze()
    A = data['A'].squeeze()
    FD_Z = data['FD_Z'].squeeze()
    R_length = data['R_length'].squeeze()
    R_width = data['R_width'].squeeze()
    R_depth = data['R_depth'].squeeze()
    FD_width = data['FD_width'].squeeze()
    FD_depth = data['FD_depth'].squeeze()
    DSC = data['DSC']
    DSU = data['DSU']

    # === STRUCTURE INIZIALIZATION ===
    OCN = {}
    OCN['thrA'] = thrA
    cellsize = 10000
    OCN['cellsize'] = cellsize
    OCN['nres'] = nReservoirs

    OCN['geometry'] = {}
    OCN['geometry']['NX'] = int(np.ceil(np.max(X / cellsize)))
    OCN['geometry']['NY'] = int(np.ceil(np.max(Y / cellsize)))
    OCN['geometry']['XI'] = np.ceil(X / cellsize).astype(int)
    OCN['geometry']['YI'] = np.ceil(Y / cellsize).astype(int)
    OCN['geometry']['SCX'] = SCX
    OCN['geometry']['SCY'] = SCY

    # SubCatchment map
    OCN['MSC'] = np.zeros((OCN['geometry']['NX'], OCN['geometry']['NY']))
    for i in range(len(X)):
        xi = OCN['geometry']['XI'][i] - 1
        yi = OCN['geometry']['YI'][i] - 1
        OCN['MSC'][xi, yi] = CTC[i]

    OCN['nNodes'] = int(np.max(CTC))
    OCN['outlet'] = outlet

    
    OCN['FD'] = {
        'A': FD_A,
        'X': FD_X,
        'Y': FD_Y,
        'downNode': FD_downNode
    }

    # Subbasin areas
    OCN['SC_AccArea'] = A
    OCN['SC_Area'] = np.zeros(OCN['nNodes'])
    for nn in range(OCN['nNodes']):
        OCN['SC_Area'][nn] = np.sum(CTC == (nn + 1)) * cellsize**2

    # Hydrological connectivity matrix
    OCN['W'] = np.zeros((OCN['nNodes'], OCN['nNodes']))
    for nn in range(OCN['nNodes']):
        temp = np.where(downNode == (nn+1))[0]
        OCN['W'][nn, temp] = 1

    # Nodes distances
    SCX_col = SCX.reshape(-1, 1)
    SCY_col = SCY.reshape(-1, 1)
    OCN['Dist'] = np.sqrt((SCX_col - SCX_col.T)**2 + (SCY_col - SCY_col.T)**2)

    # === Introducing reservoirs ===
    np.random.seed(seed_value+seedchange)
    river_pixels = np.where(FD_A >= thrA)
    nReservoirs = min(nReservoirs, len(river_pixels[0]))
    if nReservoirs > 0:
        reservoir_pixels = np.random.choice(river_pixels[0], nReservoirs, replace=False)
    else:
        reservoir_pixels = np.array([], dtype=int)

    OCN['ReservoirsPixel'] = np.zeros_like(FD_A, dtype=bool)
    OCN['ReservoirsPixel'][reservoir_pixels] = True
    OCN['ReservoirNodes'] = np.unique(CTC[reservoir_pixels]) if nReservoirs > 0 else np.array([])
    np.random.seed(seed_value)
    # === Snail habitat suitability ===
    river = (FD_A >= thrA)
    distance = np.zeros_like(FD_A)
    for i in range(FD_A.size):
        x_i, y_i = X[i], Y[i]
        dv = np.sqrt((x_i - X[river])**2 + (y_i - Y[river])**2)
        distance.flat[i] = np.min(dv)

    scaling_distance = -np.log(0.25) / 20000
    OCN['FD']['Z'] = FD_Z
    SnailHabProb = 0.8 * np.exp(-FD_Z * 0.0008) * np.exp(-scaling_distance * distance)
    
    # Effect of reservoirs
    for rp in np.where(OCN['ReservoirsPixel'].flatten())[0]:
        
        rNode = CTC[rp]
        
        upstream_nodes = np.zeros_like(FD_A, dtype=bool)
        stack = [rNode]
        while stack:
            current = stack.pop()
            
            if not upstream_nodes.flat[current - 1]:
                upstream_nodes.flat[current - 1] = True
                parents = np.where(FD_downNode == current)[0] + 1
                
                #stack.extend(parents.tolist())

        upstream_pixels = np.isin(CTC, np.where(upstream_nodes)[0] + 1)
        SnailHabProb.flat[rp] = max(SnailHabProb.flat[rp], 0.95)
        upstream_only = upstream_pixels.copy()
        upstream_only.flat[rp] = False
        SnailHabProb[upstream_only] = np.exp(-FD_Z[upstream_only] * 0.00008) * np.exp(-scaling_distance * distance[upstream_only])
    
    OCN['SnailHabitatSuitability'] = np.random.rand(*FD_Z.shape) <= SnailHabProb

    # Area snail habitat per subcatchment
    OCN['SC_SnailHabitat_Area'] = np.zeros(OCN['nNodes'])
    for nn in range(OCN['nNodes']):
        mask = (CTC == (nn + 1))
        OCN['SC_SnailHabitat_Area'][nn] = np.sum(OCN['SnailHabitatSuitability'][mask]) * cellsize**2

    # river features
    OCN['SC_RiverLength'] = R_length
    OCN['SC_RiverWidth'] = R_width
    OCN['SC_RiverDepth'] = R_depth
    OCN['SC_Volume'] = R_length * R_width * R_depth
    OCN['FD']['width'] = FD_width
    OCN['FD']['depth'] = FD_depth

    OCN['distW'] = DSC + DSC.T + DSU + DSU.T

    return OCN

