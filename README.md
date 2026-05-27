# Dams reshape water-based disease spread across river networks

This repository contains the code to reproduce all the results presented in the manuscript:

> **“Dams reshape water-based disease spread across river networks”**

By running the main file `MainScript.ipynb` cell by cell, it is possible to reproduce all the results, figures, and supplementary figures presented in the manuscript.

The *Opisthorchis viverrini* spatially explicit model is based on the framework presented in Trevisin et al. (2025). In this work, we focus on fish mobility along river corridors and, for the first time, introduce the effect of reservoirs and dam construction on *O. viverrini* spatial dynamics.

The outputs generated from the numerical solution of the differential equations are automatically saved in the `Results/` folder (mostly in `.joblib` format). These outputs are then used to reproduce the figures presented in the manuscript.

---

## Repository structure

```text
MainScript.ipynb        # Main notebook to reproduce simulations and figures
All other python files  # Python functions used throughout the project
Results/                # Simulation outputs
dataOCN/                # Synthetic river networks generated with the OCN algorithm
```

---

## Data

To run the scripts, the following elements must be available in the same working directory:

- `MainScript.ipynb`
- all required Python function files
- the folder `dataOCN/`, containing the 10 synthetic river networks generated in R using the Optimal Channel Network (OCN) algorithm

Large simulation outputs and supplementary datasets are hosted on Zenodo.

---

## Python environment

The code requires a Python environment with the required libraries needed to run the simulations and generate the figures.

Users are encouraged to create a dedicated Python environment and install the required dependencies before running the notebook.

---
