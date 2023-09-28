# MRPyCM ![WIP](https://img.shields.io/badge/-WIP-blue)
A Python implementation of the Polarizable continuum model (PCM) in [MRChem](https://github.com/MRChemSoft/mrchem).
This is implemented using [VAMPYR](https://github.com/MRChemSoft/vampyr).
This is a work in progress.


## Installation
### Building from source

[!IMPORTANT]
(python 3.11.4+ and pip 23.2.1+ required)

You can install add the package to your python/pip/conda environment by cloning the repository
```bash
git clone git@github.com:Gabrielgerez/MRPyCM.git
```
and then installing it with pip 

```bash
cd MRPyCM
pip install .
``` 

## Usage
The code can be used as an import library by:
``` python
import MRPycm as mpcm
```
or you can run the executable `computeSolvent.py` as:

```bash
python computeSolvent.py <dict_str>
```

where `dict_str` is a dictionary of input parameters as

```
{
    'order': int,                           # Polynomial order.
    'box': list of floats,                  # Simulation box size.
    'prec': float,                          # Apply and convergence precision.
    'cav_coords': list of list of floats,   # Coordinate of the center of each sphere.
    'radii': list of floats,                # Radii of each sphere.
    'sigma': float,                         # Width of the boundary of each sphere.
    'eps_out': float,                       # Value of the relative permittivity of the solvent at T=298.15K.
    'perm_type': str,                       # Formulation of the permittivity either "exponential" or "linear".
    'solvent_type: str,                     # Which solver to use, avaliable are "standard" (GPESolver), "pb" (PBSolver) and "lpb" (LPBSolver).
    'I': float                              # Ionic strength in mol/L. Can also be interpreted as the electrolyte concentration.
    'max_iter': int,                        # Maximum number of iterations before exiting the solver.
    'kain_hist': int,                       # Number of iterations to keep in the history of the Kain solver.
}
```
