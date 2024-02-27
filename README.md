# Activation Volume Optimizer

This repo contains a set of python scripts, shell scripts, and a version of the OPLD code to determine the activation energies of different sized clusters in a pure crystalline atomic structure.

## Methodology for Obtaining Activation Energies

To obtain the average activation energies of a cluster or a point defect three steps must be taken.

1. Create multiple configurations of the defect
  Use the 'configuration_finder' script to accomplish this. Use -h as an argument to see command line options.

2. Run dimer searches of the configurations to find possible activation energies.
  The OPLD executable must be compiled first. See How To for OPLD dimer searches below.

3. Create the histogram from the dimer search results.
  Use the 'plotting' script to create the histogram. Use -h to see required arguments.

**Use the 'run' bash shell script to automate the processes of obtaining the activation energy**
Before using the bash script you must first compile the OPLD executable.

## What You Will Need
You just need one thing, a LAMMPS script which builds a perfect crystal. Make sure that the crystal is equilibrated through the script.

## Compiling OPLD executable

To compile OPLD, you must first build [LAMMPS](https://www.lammps.org/#gsc.tab=0) version stable_29Sep2021. This can be done by cloning LAMMPS from Git and checking out the appropiate tag. You want to build the static library with all the necessary LAMMPS packages.
Don't forget to install the LAMMPS Python package as well.

Head in to the 'opldSrc' directory. Run the make file by simply calling the Makefile. Be sure to compile OPLD with the same OpenMPI version used in the LAMMPS compilation process.

## How To OPLD
The OPLD executable requires two arguments, the path to a LAMMPS script which creates a perfect crystal and the path a configuration file created by the 'configuration_finder' script. Activation energies will be printing to stdout(screen).