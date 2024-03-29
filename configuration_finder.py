# For the requested defect type, this script will calculate activation
# energies, energy per atom, and hydrostatic stress as a function of
# AV radius. It will repeat this for each configuration of the defect
# type. The output are the following graphs:
#  - Activation energy vs AV radius (each configuration)
#  - Energy per atom vs AV radius (each configuration)
#  - Hydrostatic stress vs AV radius (each configuration)
#  - Activation Energy histogram (each configuration)

import argparse
from lammps import lammps
from mpi4py import MPI
from avo import Avo

# Initialize MPI
comm = MPI.COMM_WORLD

# print("MPI version:", MPI.Get_library_version())

# parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("defect_size", type=int, help="size of defect to analyze")
parser.add_argument(
    "--output_dir",
    type=str,
    help="output directory otherwise ./tempData will be created",
)
parser.add_argument(
    "--data_dir",
    type=str,
    help="data directory ready for saddle point searching (skips condfiguration generation)",
)
parser.add_argument(
    "lammps_filename", type=str, help="lammps input file for perfect crystal"
)
parser.add_argument(
    "lattice_parameter", type=float, help="lattice parameter of the crystal"
)
parser.add_argument(
    "--number_configurations",
    type=int,
    default=5,
    help="number of configurations to analyze",
)
parser.add_argument(
    "--equilibration_temp",
    type=float,
    default=800.0,
    help="temperature (K) to equilibrate to",
)

args = parser.parse_args()

defect_size = args.defect_size
output_dir = args.output_dir
data_dir = args.data_dir
lammps_filename = args.lammps_filename
lattice_parameter = args.lattice_parameter
number_configs = args.number_configurations
equilibration_temp = args.equilibration_temp

avo = Avo(
    mpiComm=comm,
    resultsDir=output_dir,
    dataDir=data_dir,
)
# First determine the best possible configurations of the defect type
avo.findConfigurations(
    defect_size,
    number_configs,
    lammps_filename,
    lattice_parameter,
    equilibration_temp,
)

if avo.rank == 0:
    avo.createOPLDInputFile()
