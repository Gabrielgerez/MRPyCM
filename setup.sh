#!/bin/bash

# update external libraries
echo "Updating external libraries"
git submodule update --remote --recursive --init

# Install dependencies
# conda is a requirement for installing all dependencies
eval "$(conda shell.bash hook)"

find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

# check if conda environment exists, if not create it.
if find_in_conda_env ".*MRPyCM-env.*" ; then
    echo "MRPyCM-env found"
else 
    echo "MRPyCM-env not found, creating it"
    conda env create -f environment.yml 
fi

# activate conda environment
echo "Activating MRPyCM-env"
conda activate MRPyCM-env

# install MRPyCM
echo "Installing MRPyCM"
pip install .

# run tests
echo "Running tests"
pytest --ignore=external