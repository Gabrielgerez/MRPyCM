# MRPyCM ![WIP](https://img.shields.io/badge/-WIP-blue)
A Python implementation of the Polarizable continuum model (PCM) in[MRChem](https://github.com/MRChemSoft/mrchem).
This is implemented using [VAMPYR](https://github.com/MRChemSoft/vampyr).
This is a work in progress.


## installation
To use the code you download it from github and run the script `setup.sh` as follows:
```bash
git clone git@github.com:Gabrielgerez/MRPyCM.git
cd MRPyCM
./setup.sh
``` 

## Usage
The code can be used as an import library by:
``` python
import MRPycm as mpcm
```
or you can run it as an executable through the command line by:

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
}
```
