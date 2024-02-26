#!/bin/bash

# Author: Angel Chavira
# Date: 2/23/2024

# This script is used to run the dimer searches for all the configurations found in the
# given directory ending with the '.lammpstrj' suffix. 

# The configuration _finder.py python script should be first be called to create all the 
# configurations wanted to be used in the dimer searches.

script_path=$(realpath $0)
opld_executable=$(dirname $script_path)/opldSrc/opld

# Ask where to find certain critical files
echo "Which directory would you like to run the dimer searches in?"
read directory
if [ ! -d "${directory}" ]; then
    echo "Directory ${directory} does not exist"
    exit 1
else
    echo "Using configurations found in ${directory}"
fi

echo "Where is the lammps input file located?"
read lammps_input
if [ ! -f "${lammps_input}" ]; then
    echo "File ${lammps_input} does not exist"
    exit 1
else
    echo "Using lammps input file located at ${lammps_input}"
fi

handle_interrupt() {
    echo -e "\nInterrupt received, exiting..."
    exit 1
}
trap handle_interrupt SIGINT

# Do the dimer searches
for file in ${directory}/*.lammpstrj
do
    echo "Running dimer search for ${file}" >> dimer_results.log
    mpirun ${opld_executable} ${lammps_input} ${file} >> ${directory}/dimer_results.log
done