import numpy as np

def zipf(rank, expn=1, minP=250):
    """
    Generate a Zipf-distributed array.
    
    rank : array-like
        The ranks of the elements (integers or array).
    expn : float
        The exponent of the Zipf distribution.
    minP : float
        Scaling factor (minimum value after normalization).
        
    Returns
    -------
    P : np.ndarray
        Zipf-distributed vector.
    """
    rank = np.array(rank+1, dtype=float)
    
    H = np.sum(1.0 / (rank ** expn))
    print(H)
    P = 1.0 / (rank ** expn) / H   # normalized probabilities
    P = P / np.min(P) * minP       # scale so min(P) = minP
    return P
