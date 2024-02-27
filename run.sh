#!/bin/bash

# Author: Angel Chavira
# Date: 2/23/2024

# This script is used to run the entire activation optimization processes without
# manunally running each step in the process. Use this script if you are starting
# from scratch. This will only find the activation energy for a single cluster size.

# ************ DO NOT MOVE THIS FILE; IT COULD BREAK CODEBASE ************

# Find where the current script is located. Need this to find the OPLD executable.

script_path=$(realpath $0)
script_dir=$(dirname $script_path)
dimer_search_script=${script_dir}/dimer_searches.sh
config_finder=${script_dir}/configuration_finder.py
histo_plotter=${script_dir}/plotting.py

# Ask critical questions
echo "What sized cluster would you like to find the activation energy for?"
read cluster_size
echo "Using cluster size ${cluster_size}"

echo "Where is the lammps input file located?"
read lammps_input

echo "What is the lattice constant?"
read lattice_constant

echo "What is the temperature?"
read temperature

configurations_dir=./configurations${cluster_size}

handle_interrupt() {
    echo -e "\nInterrupt received, exiting..."
    exit 1
}
trap handle_interrupt SIGINT

# Find various configurations
mpirun python3 ${config_finder} ${cluster_size} ${lammps_input} ${lattice_constant} --output_dir ${configurations_dir} \
    --number_configurations 5 --equilibration_temp ${temperature}

# Do dimer searches for all configurations in the configurations directory
echo -e "${configurations_dir}\n${lammps_input}" | ${dimer_search_script}

# Plot the activation energies
python3 ${histo_plotter} ${configurations_dir} ${temperature}