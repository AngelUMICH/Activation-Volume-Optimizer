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
_randomDisplacement = 0.2 # lattice parameter
_lammpsTimestep = 0.001  # in picoseconds (10^-15)
_lammpsTimePerSnapshot = 3  # in picoseconds (10^-15)
_lammpsTimestepsPerSnapshot = _lammpsTimePerSnapshot / _lammpsTimestep


# ********************************************************************#
class _MDRunner:

    def __init__(self, mpiComm: MPI.Comm):
        super().__init__()
        self.lammps = lammps(cmdargs=["-log", "none"])
        self.mpiComm = mpiComm
        self.rank = self.mpiComm.Get_rank()
        self.mpiSize = self.mpiComm.Get_size()

    def __setupLammps(self):
        # Needed to get scatter_atoms to work
        self.lammps.command("atom_modify map array")
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
        self.lammps.command(f"dump referenceDump all xyz 1 ./pos_ref.dat")
        self.lammps.command("run 0")
        self.lammps.command("undump referenceDump")

    def __createDefect(self):
        # Either delete or add atoms to create the defect
        if self.rank == 0:
            print(f"Creating defect of size {self.defectSize}")

        if self.defectSize > 0:
            self.__addAtoms(self.defectSize)
            if self.rank == 0:
                print(f"Created interstitial defect of size {self.defectSize}")
        elif self.defectSize < 0:
            self.__deleteAtoms(abs(self.defectSize))
            if self.rank == 0:
                print(f"Deleted vacancy defect of size {-self.defectSize}")

    def __addAtoms(self, numberAtomsAdd: int):
        # Add atoms to create a defect cluster
        nlocal = self.lammps.extract_global("nlocal") # used to select random atom type
        atomPositions = self.lammps.gather_atoms("x", 1, 3)
        closestAtomIds = self.__findClosestAtoms(atomPositions, numberAtomsAdd)
        # Pick a random <100> dumbbell direction for ints in Nickel
        orientation = [0, 0, 0]
        orientation[random.randint(0, 2)] = 1
        orientation = self.mpiComm.bcast(orientation, root=0) # synchronize orientation
        if self.rank == 0:
            print(f"Orientation of dumbbell: {orientation}")
            print("Modify code for orientation other than <100>")
        # New atom positions
        newAtoms = []
        for id in closestAtomIds:
            randomAtomId = random.randint(0, nlocal - 1)
            randomType = self.atomTypes[randomAtomId]
            # Find location of new atom
            x = atomPositions[id * 3 + 0] + _randomDisplacement * self.latticeParameter * orientation[0]
            y = atomPositions[id * 3 + 1] + _randomDisplacement * self.latticeParameter * orientation[1]
            z = atomPositions[id * 3 + 2] + _randomDisplacement * self.latticeParameter * orientation[2]
            newAtoms.append((randomType, x, y, z))
            # Find new location of original atom
            atomPositions[id * 3 + 0] = x - _randomDisplacement * self.latticeParameter *  orientation[0] * 2
            atomPositions[id * 3 + 1] = y - _randomDisplacement * self.latticeParameter *  orientation[1] * 2
            atomPositions[id * 3 + 2] = z - _randomDisplacement * self.latticeParameter *  orientation[2] * 2

        # Scatter new atom positions
        self.lammps.scatter_atoms("x", 1, 3, atomPositions)
        # Create the new atom. Create_atoms must be after scatter_atoms
        for atom in newAtoms:
            randomType, x, y, z = atom
            self.lammps.command(
                "create_atoms %d single %f %f %f units box" % (randomType, x, y, z)
            )

    def __deleteAtoms(self, numberAtomDelete: int):
        # Update the MPI ranks with correct data
        self.crystalCenterPosition = self.mpiComm.bcast(
            self.crystalCenterPosition, root=0
        )
        # Delete atoms which are near crystal center
        self._createDeleteGroup(numberAtomDelete)
        self.lammps.command(f"delete_atoms group deleteGroup")
        self.lammps.command("group deleteGroup clear")

    def _createDeleteGroup(self, numberAtomDelete: int):
        # Create a group of atoms within the delete sphere region
        deleteGroupIds = []
        atomPositions = self.lammps.gather_atoms("x", 1, 3)

        if self.rank == 0:
            deleteGroupIds = self.__findClosestAtoms(atomPositions, numberAtomDelete)

        # Bcast ids and create group
        deleteGroupIds = self.mpiComm.bcast(deleteGroupIds, root=0)
        for id in deleteGroupIds:
            self.lammps.command(f"group deleteGroup id {id}")

    def __findClosestAtoms(self, positions, natoms) -> list:
        # Find the closest natoms to the crystal center
        closestAtoms = [
            (0, 1_000_000) for _ in range(natoms)
        ]  # (a,b) a= id, b= distance
        natoms = self.lammps.get_natoms()

        for atomId in range(natoms):
            x = positions[atomId * 3 + 0] / self.latticeParameter
            y = positions[atomId * 3 + 1] / self.latticeParameter
            z = positions[atomId * 3 + 2] / self.latticeParameter
            deltaRR = (
                (x - self.crystalCenterPosition[0])
                * (x - self.crystalCenterPosition[0])
                + (y - self.crystalCenterPosition[1])
                * (y - self.crystalCenterPosition[1])
                + (z - self.crystalCenterPosition[2])
                * (z - self.crystalCenterPosition[2])
            )
            if deltaRR < closestAtoms[-1][1]:
                closestAtoms[-1] = (atomId + 1, deltaRR)
                closestAtoms.sort(key=lambda x: x[1])

        # Return a list of only the ids of closestAtoms list
        return [atom[0] for atom in closestAtoms]

    def __runDynamics(self):
        if self.rank == 0:
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
        self.lammps.command("reset_timestep 1")
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
        self.numberConfigs = numberConfigs
        self.inputLammpsFilename = inputLammpsScript
        self.latticeParameter = latticeParameter
        self.equilibrationTemp = equilibrationTemp
        self.__setupLammps()
        self.__createReferenceConfiguration()
        if self.configurationState:
            self.__createDefect()
            self.__runDynamics()
        self.lammps.close()
