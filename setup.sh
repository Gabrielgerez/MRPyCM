#!/bin/bash

# Install dependencies
# conda is a requirement for installing all dependencies
find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

eval "$(conda shell.bash hook)"

# create or activate conda environment
## check if environment exists, if not create it then activate it
if find_in_conda_env ".*MRPyCM-env.*" ; then
    echo "MRPyCM-env found, activating it"
    conda activate MRPyCM-env
else 
    echo "MRPyCM-env not found, creating it"
    conda env create -f environment.yml 
    echo "Activating MRPyCM-env"
    conda activate MRPyCM-env
fi

# update external libraries
echo "Updating external libraries"
git submodule update --remote --recursive --init

# # build vampyr
echo "Building vampyr"
cd external/vampyr
pip install .