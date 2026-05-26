import numpy as np
import math
from zipf import zipf
from calculateWRecursive import calculateWRecursive

def build_setup_Reservoirs(OCN, p, TotalPopulation, 
                           DownstreamAccumulation=False, 
                           Unify=False, 
                           seed=None):
    """
    Python translation of MATLAB build_setup_Reservoirs.
    Reproduces identical logic and output structure.
    """

    if seed is not None and not np.isnan(seed):
        np.random.seed(int(seed))

    s = {}
    s['nNodes'] = OCN['nNodes']
    s['par'] = p

    # Hydrological connectivity
    s['X'] = OCN['W'].copy()

    out = np.argmax(OCN['SC_AccArea'])
    if OCN['SC_RiverLength'][out] == 0:
        OCN['SC_RiverLength'][out] = OCN['cellsize']

    # === POPULATIONS ===
    # Generate random human population Zipf distributed
    if DownstreamAccumulation:
        rankacc = np.argsort(OCN['SC_AccArea'])[::-1]
        ID = np.zeros_like(rankacc)
        for rank in range(OCN['nNodes']):
            ID[rankacc[rank]] = rank + 1
        s['H'] = zipf(ID, 1, 250)
    else:
        s['H'] = zipf(np.random.permutation(OCN['nNodes']), 1, 250)

    s['H'] = s['H'] / np.sum(s['H']) * TotalPopulation

    # Fish population
    s['F'] = p['dF'] * (OCN['SC_RiverLength'] * OCN['SC_RiverWidth'])

    # Habitat area
    s['A'] = np.copy(OCN['SC_SnailHabitat_Area'])
    s['A'][s['A'] == 0] = np.nan

    # Snail population
    s['S'] = OCN['SC_SnailHabitat_Area'] * p['dS']

    # Fish mobility matrix
    visited = np.zeros(s['nNodes'], dtype=bool)
    s['W'] = calculateWRecursive(0, s['X'] + s['X'].T, s['F'], np.zeros_like(s['X']), visited, 1)
    print(np.mean(s['W'][s['W'] > 0]))
    if math.isnan(np.mean(s['W'][s['W'] > 0])):
        print('why are we here')
    else:
        s['W'] = s['W'] * s['par']['lambda_F'] / np.mean(s['W'][s['W'] > 0])
    

    # === FISH MARKET === Not included in our study
    s['chi'] = s['par']['c'] * s['H']
    s['par']['U'] = np.nansum(s['chi'] * s['F']) / np.nansum(s['H'])
    s['delta'] = -s['chi'] * s['F'] + s['par']['U'] * s['H']
    s['delta'][s['delta'] < 0] = 0
    s['sigma'] = s['chi'] * s['F'] - s['par']['U'] * s['H']
    s['sigma'][s['sigma'] < 0] = 0

    # Snail exposure reduction
    s['par']['xi'] = np.random.rand(len(s['H'])) * np.exp(
        -np.random.rand(len(s['H'])) * 0.00001 * s['H']
    )

    # Raw fish quota
    s['par']['epsilon'] = s['par']['epsilon0'] * np.ones(s['nNodes'])

    # Fish exposure
    s['par']['theta'] = (
        s['par']['rho_C'] * s['par']['beta_C'] * s['S'] / s['A'] /
        (s['par']['mu_C'] + s['par']['beta_C'] * s['F'] / s['A'])
    )
    s['par']['theta'][np.isnan(s['par']['theta'])] = 0

    # === UNIFY ===
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
        s['par']['epsilon'] = np.nansum(s['par']['epsilon'] * s['H']) / np.nansum(s['H'])
        s['H'] = np.nansum(s['H'])
        s['S'] = np.nansum(s['S'])
        s['F'] = np.nansum(s['F'])
        s['par']['theta'] = (
            s['par']['rho_C'] * s['par']['beta_C'] * s['S'] / s['A'] /
            (s['par']['mu_C'] + s['par']['beta_C'] * s['F'] / s['A'])
        )

    return s
