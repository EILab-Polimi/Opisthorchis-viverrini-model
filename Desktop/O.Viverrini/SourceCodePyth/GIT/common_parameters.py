import numpy as np

def common_parameters():
    """
    Define and return a dictionary of common parameters.
    This is a direct translation of the MATLAB function common_parameters.m
    """

    par = {}

    # Mortality
    par["mu_H"] = 1 / (365.25 * 68.5)       # Life expectancy (Laos)
    par["mu_W"] = 1 / (365.25 * 10)         # Worm mortality
    par["mu_E"] = -np.log(0.8) / 21         # Egg mortality

    # Egg shedding
    par["rho_E"] = 4.848252470797727e+02 * 200   # Crellen (multiplied by 200 grams of stool)
    par["alpha"] = 23.241262353988635            # Crellen (10 for 1 worm infection)

    # Cercariae
    par["rho_C"] = (263 * 6 + 98 * 9 + 91 * 17) / (6 + 9 + 17)  # cercariae/day
    par["mu_C"] = -np.log(0.01) / 2

    # Common exposure to raw fish
    par["epsilon0"] = 0.2

    # Calibrated parameters
    par["beta_E"] = 5e-9
    par["beta_C"] = 0.04     # originally 0.005 in MATLAB
    par["c"]      = 5e-9
    par["dF"]     = 0.08
    par["dS"]     = 20       # originally 10^(3.3820) in MATLAB
    par["mu_F"]   = 1 / (365 * 2.5)
    par["mu_S"]   = 1 / (365 * 1)

    # Mobility parameters
    par["D"] = 5e5
    par["lambda_F"] = 0.01

    return par
