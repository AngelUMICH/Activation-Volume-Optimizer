from lammps import lammps
import random
from mpi4py import MPI

# Variables needed for the tools module
_xhi = 3
_xlo = 0
_yhi = 4
_ylo = 1
_zhi = 5
_zlo = 2
_radiusBuffer = 0.05
_lammpsTimestep = 0.001  # in picoseconds (10^-15)
_lammpsTimePerSnapshot = 3  # in picoseconds (10^-15)
_lammpsTimestepsPerSnapshot = _lammpsTimePerSnapshot / _lammpsTimestep
# ********************************************************************#
class _MDRunner:

    def __init__(self, mpiComm: MPI.Comm):
        super().__init__()
        self.lammps = lammps(cmdargs=["-screen", "none", "-log", "none"])
        self.mpiComm = mpiComm
        self.rank = self.mpiComm.Get_rank()
        self.mpiSize = self.mpiComm.Get_size()

    def __setupLammps(self):
        # inputLammpsFilename is the filename of a minimum working lammps script
        self.lammps.file(self.inputLammpsFilename)
        self.atomTypes = self.lammps.extract_atom("type")
        self.crystalNumberAtoms = self.lammps.get_natoms()
        self.crystalBoxBounds = self.lammps.extract_box()
        self.crystalBoxBounds = [
            item for sublist in self.crystalBoxBounds[0:2] for item in sublist
        ]
        xCenter = (
            self.crystalBoxBounds[_xlo]
            + (self.crystalBoxBounds[_xhi] - self.crystalBoxBounds[_xlo])
            / 2
            / self.latticeParameter
        )
        yCenter = (
            self.crystalBoxBounds[_ylo]
            + (self.crystalBoxBounds[_yhi] - self.crystalBoxBounds[_ylo])
            / 2
            / self.latticeParameter
        )
        zCenter = (
            self.crystalBoxBounds[_zlo]
            + (self.crystalBoxBounds[_zhi] - self.crystalBoxBounds[_zlo])
            / 2
            / self.latticeParameter
        )
        self.crystalCenterPosition = (xCenter, yCenter, zCenter)

        # setup basic lammps commands
        self.lammps.command(f"timestep %f" % (_lammpsTimestep))
        self.totalNumberTimesteps = self.numberConfigs * _lammpsTimestepsPerSnapshot

    def __createReferenceConfiguration(self):
        # Create a reference of the perfect crystal by dumping
        # the configuration at the beginning of the simulation
        self.lammps.command(
            f"dump referenceDump all xyz 1 ./pos_ref.dat")
        self.lammps.command("run 0")
        self.lammps.command("undump referenceDump")

    def __createDefect(self):
        # Either delete or add atoms to create the defect
        print(f"Creating defect of size {self.defectSize}")
        if self.defectSize > 0:
            self.__addAtoms(self.defectSize)
            print(f"Created interstitial defect of size {self.defectSize}")
        elif self.defectSize < 0:
            self.__deleteAtoms(self.defectSize)
            print(f"Deleted vacancy defect of size {-self.defectSize}")

    def __addAtoms(self, numberAtomsAdd: int):
        # Add atoms to create a defect cluster
        for _ in range(numberAtomsAdd):
            if self.rank == 0:
                nlocal = self.lammps.extract_global("nlocal")
                randomAtomId = random.randint(0, nlocal - 1)
                randomType = self.atomTypes[randomAtomId]
                x = (
                    random.random() - 0.5
                ) * self.defectRegionRadius * 2 + self.crystalCenterPosition[0]
                y = (
                    random.random() - 0.5
                ) * self.defectRegionRadius * 2 + self.crystalCenterPosition[1]
                z = (
                    random.random() - 0.5
                ) * self.defectRegionRadius * 2 + self.crystalCenterPosition[2]
            elif self.rank != 0:
                randomType = None
                x = None
                y = None
                z = None
            randomType = self.mpiComm.bcast(randomType, root=0)
            x = self.mpiComm.bcast(x, root=0)
            y = self.mpiComm.bcast(y, root=0)
            z = self.mpiComm.bcast(z, root=0)
            print(
                f"Rank {self.rank} selected type {randomType} at position {x} {y} {z}"
            )
            self.lammps.command(
                "create_atoms %d single %f %f %f units lattice" % (randomType, x, y, z)
            )

    def __deleteAtoms(self, numberAtomDelete: int):
        # Delete atoms which are near crystal center
        self.lammps.command(
            f"region defectRegion sphere %f %f %f %f units lattice"
            % (
                self.crystalCenterPosition[0],
                self.crystalCenterPosition[1],
                self.crystalCenterPosition[2],
                self.defectRegionRadius,
            )
        )
        self.lammps.command("group defectGroup region defectRegion")
        self.lammps.command(
            f"delete_atoms random count %d yes defectGroup defectRegion %d"
            % (numberAtomDelete, random.randint(1, 1000000))
        )

    def __runDynamics(self):
        print(
            f"Running dynamics for {self.totalNumberTimesteps} timesteps to find %d configurations"
            % (self.numberConfigs)
        )
        # First minimize the initial configuration
        self.lammps.command("minimize 1.0e-8 1.0e-8 1000 10000")
        # Equilibrate the system to the desired temperature
        self.lammps.command(
            f"velocity all create %f %d"
            % (self.equilibrationTemp, random.randint(1, 1000000))
        )
        self.lammps.command(
            "fix 1 all nvt temp %f %f 0.1"
            % (self.equilibrationTemp, self.equilibrationTemp)
        )
        self.lammps.command("run 2000")

        # Now setup the dumpfiles of each snapshot for dimer searches
        self.lammps.command("reset_timestep 0")
        self.lammps.command(
            f"dump snapshotDump all xyz %d %s/configuration*.lammpstrj"
            % (_lammpsTimestepsPerSnapshot, self.dataFolder)
        )

        # Run the dynamics for the desired number of timesteps
        self.lammps.command(f"run %d" % (self.totalNumberTimesteps))
        self.lammps.command("undump snapshotDump")
        self.lammps.command("unfix 1")

    def findConfigurations(
        self,
        defectSize: int,
        numberConfigs: int,
        inputLammpsScript: str,
        latticeParameter: float,
        equilibrationTemp: float,
    ):
        # This function will find the best configurations of the defect type
        self.defectSize = defectSize
        self.defectRegionRadius = self.defectSize ** (1/3) + _radiusBuffer
        self.numberConfigs = numberConfigs
        self.inputLammpsFilename = inputLammpsScript
        self.latticeParameter = latticeParameter
        self.equilibrationTemp = equilibrationTemp
        self.__setupLammps()
        self.__createReferenceConfiguration()
        if (self.configurationState):
            self.__createDefect()
            self.__runDynamics()
        self.lammps.close()

# ********************************************************************#