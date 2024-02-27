import os
from mpi4py import MPI
from mdrunner import _MDRunner
from opldrunner import _OPLDRunner


# ********************************************************************#
# This class will house all the necessary tools for the activation relaxation optimizer
class Avo(_MDRunner, _OPLDRunner):

    def __init__(self, mpiComm: MPI.Comm, resultsDir: str, dataDir: str):
        super().__init__(mpiComm)
        self.__determineAvoState(dataDir)
        self.__CreateDataFolder(resultsDir, dataDir)

    def __determineAvoState(self, dataDir: str):
        if dataDir:
            self.configurationState = False
        elif not dataDir:
            self.configurationState = True

    def __CreateDataFolder(self, outputDir: str, inputDir: str):
        # Check to see if we even need to create a data folder
        if inputDir:
            self.dataFolder = inputDir
        else:
            # Create a temporary folder to store the data
            if not outputDir:
                self.dataFolder = os.getcwd() + "/tempData"
            else:
                self.dataFolder = os.path.abspath(outputDir)
            if not os.path.exists(self.dataFolder):
                if self.mpiComm.Get_rank() == 0:
                    print(f"Creating data folder {self.dataFolder}")
                    os.makedirs(self.dataFolder)
