import numpy as np
def calculateWRecursive(node, H, F, W, visited, val):
    connectedNodes = np.where(H[node, :] == 1)[0]
    W[connectedNodes,node] = val
    visited[node] = True
    for i in connectedNodes:
        #W[i, node] = val
        if not visited[i]:
            W[node, i] = W[i, node] * F[node] / F[i]
            val = W[i, node] * F[node] / F[i]
            W = calculateWRecursive(i, H, F, W, visited, val)
            visited[i] = True
    return W
