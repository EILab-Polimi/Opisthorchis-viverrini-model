import numpy as np
import math
from zipf import zipf
from calculateWRecursive import calculateWRecursive
def build_setup(OCN, p, TotalPopulation, DownstreamAccumulation=False, Unify=False, seed=None):
    """
    Build simulation setup from OCN structure and parameters.
    """

    if seed is not None:
        np.random.seed(seed)

    s = {}
    s['nNodes'] = OCN['nNodes']
    s['par'] = p.copy()

    # Hydrological connectivity
    s['X'] = OCN['W'].copy()
    out = np.argmax(OCN['SC_AccArea'])
    if OCN['SC_RiverLength'][out] == 0:
        OCN['SC_RiverLength'][out] = OCN['cellsize']

    # --- Human population ---
    if DownstreamAccumulation:
        rankacc = np.argsort(-OCN['SC_AccArea'])
        ID = np.zeros(OCN['nNodes'])
        for rank, node in enumerate(rankacc):
            ID[node] = rank + 1
        s['H'] = zipf(ID, 1, 250)
    else:
        s['H'] = zipf(np.random.permutation(OCN['nNodes']), 1, 250)
    
    s['H'] = s['H'] / np.sum(s['H']) * TotalPopulation

    # --- Fish population ---
    s['F'] = p['dF'] * (OCN['SC_RiverLength'] * OCN['SC_RiverWidth'])

    # --- Snail population ---
    s['A'] = OCN['SC_SnailHabitat_Area'].copy()
    s['A'][s['A'] == 0] = np.nan
    s['S'] = OCN['SC_SnailHabitat_Area'] * p['dS']

    # --- Compute fish mobility matrix ---
    visited = np.zeros(s['nNodes'], dtype=bool)
    s['W'] = calculateWRecursive(0, s['X'] + s['X'].T, s['F'], np.zeros_like(s['X']), visited, 1)
    print(np.mean(s['W'][s['W'] > 0]))
    if math.isnan(np.mean(s['W'][s['W'] > 0])):
        print('why are we here')
    else:
        s['W'] = s['W'] * s['par']['lambda_F'] / np.mean(s['W'][s['W'] > 0])

    # --- Fish market --- not included in our study
    s['chi'] = s['par']['c'] * s['H']
    s['par']['U'] = np.nansum(s['chi'] * s['F']) / np.nansum(s['H'])
    s['delta'] = -s['chi'] * s['F'] + s['par']['U'] * s['H']
    s['delta'][s['delta'] < 0] = 0
    s['sigma'] = s['chi'] * s['F'] - s['par']['U'] * s['H']
    s['sigma'][s['sigma'] < 0] = 0

    

    # Snail exposure reduction
    s['par']['xi'] = np.random.rand(s['H'].size) * np.exp(-np.random.rand(s['H'].size) * 0.00001 * s['H'])
    s['par']['epsilon'] = np.full(s['nNodes'], s['par']['epsilon0'])
    s['par']['theta'] = p['rho_C'] * p['beta_C'] * s['S'] / s['A'] / (p['mu_C'] + p['beta_C'] * s['F'] / s['A'])
    s['par']['theta'] = np.nan_to_num(s['par']['theta'])

    # --- Unify option ---
    if Unify:
        s['chi'] = s['par']['c'] * np.nansum(s['H'] * s['F']) / np.nansum(s['F'])
        s['delta'] = s['chi'] * np.nansum(s['F']) - s['par']['U'] * np.nansum(s['H'])
        s['delta'] = np.maximum(s['delta'], 0)
        s['sigma'] = -s['chi'] * np.nansum(s['F']) + s['par']['U'] * np.nansum(s['H'])
        s['sigma'] = np.maximum(s['sigma'], 0)
        s['T'] = min(s['chi'] * np.nansum(s['F']), s['par']['U'] * np.nansum(s['H']))

        s['nNodes'] = 1
        s['W'] = 0
        s['A'] = np.nansum(s['A'])
        s['par']['xi'] = np.nansum(s['par']['xi'] * s['H']) / np.nansum(s['H'])
        #s['par']['epsilon'] = np.nansum(s['par']['epsilon'] * s['H']) / np.nansum(s['H'])
        s['H'] = np.nansum(s['H'])
        s['S'] = np.nansum(s['S'])
        s['F'] = np.nansum(s['F'])
        s['par']['theta'] = p['rho_C'] * p['beta_C'] * s['S'] / s['A'] / (p['mu_C'] + p['beta_C'] * s['F'] / s['A'])

    return s
