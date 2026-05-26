import numpy as np
import networkx as nx

def build_connectivity_matrix(river, nR, dam_pass, alpha, pixel2rivnode_full, nres):
    """
    Python equivalent of the MATLAB function build_connectivity_matrix.
    """
    river = np.array(river)
    np.random.seed(nres)
    # Tutti i pixel unici
    allPixels = np.unique(river[:, 0])
    pixel2rivnode = np.array(pixel2rivnode_full)[allPixels.astype(int)]
    
    nNodes = len(allPixels)
    '''

   # Mapping pixelID → consecutive index (0-based)
    all_edges = np.concatenate((river[:, 0], river[:, 1]))
    _, _, newIdx = np.unique(all_edges, return_index=True, return_inverse=True)

    s = newIdx[:len(river)] 
        # sources
    t = newIdx[len(river):] 
    print(t)     # targets
    w = river[:, 2]              # weights (lengths)
    reservoirs = river[:, 3]     # dam flag 0/1

    # Graph construction
    G = nx.DiGraph()
    for i in range(len(river)):
        G.add_edge(int(s[i]), int(t[i]), weight=float(w[i]))
    '''
    all_edges = np.concatenate((river[:, 0], river[:, 1]))
    _, _, newIdx = np.unique(all_edges, return_index=True, return_inverse=True)

    s = newIdx[:len(river)]
    t = newIdx[len(river):]
    w = river[:, 2]
    reservoirs = river[:, 3]

    # Building the graph as in Matlab
    G = nx.Graph()  # not oriented
    G.add_weighted_edges_from(zip(s, t, w))
    unique_pixels = np.unique(all_edges)
    reservoirs_pixel = np.zeros(len(unique_pixels))
    for idx, pix in enumerate(unique_pixels):
        if np.any(river[(river[:, 0] == pix) | (river[:, 1] == pix), 3] == 1):
            reservoirs_pixel[idx] = 1
    LambdaF_mat = np.zeros((nNodes, nNodes))

    # Computing LambdaF_mat (Connectivity matrix pixel→pixel)
    for i in range(nNodes):
        i = i+1
        for j in range(nNodes):
            j=j+1
            if i == j:
                LambdaF_mat[i-1, j-1] = 1.0
            else:
                path = nx.shortest_path(G, source=i, target=j, weight='weight')
                
                dist_ij = nx.shortest_path_length(G, source=i, target=j, weight='weight')

                
            

                # Passability along the path 
                pass_ij = 1
                #for k in path[:-1]:  
                for k in path:
                    if reservoirs_pixel[k] == 1:
                        if dam_pass < 0.26:
                            #print('dam passability correct and equal 0')
                            pass_ij *= dam_pass
                        elif dam_pass == 1:
                            #print('dam passability wrong and equal 1')
                            pass_ij *= dam_pass
                        else:
                            # Generate a random passability value between 0.1 and 0.5
                            #print('dam passability wrong and equal 0.3')
                            dam_pass = np.random.uniform(0.1, 0.5)
                            pass_ij *= dam_pass
                        

                w_ij = np.exp(-alpha * dist_ij)
                LambdaF_mat[i-1, j-1] = w_ij * pass_ij

    # Computing LambdaRiver (aggreagted on river nodes)
    LambdaRiver = np.zeros((nR, nR))

    for a in range(nR):
        a=a+1
        apix = np.where(pixel2rivnode == a)[0]
        
        for b in range(nR):
            b=b+1
            bpix = np.where(pixel2rivnode == b)[0]

            if len(apix) == 0 or len(bpix) == 0:
                LambdaRiver[a-1, b-1] = 0
                continue

            sub = LambdaF_mat[np.ix_(apix, bpix)]
            if np.min(sub) == 0:
                LambdaRiver[a-1, b-1] = 0
            else:
                LambdaRiver[a-1, b-1] = np.mean(sub)

    return LambdaRiver, LambdaF_mat



