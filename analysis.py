import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

BOLTZMANN_CONSTANT = 8.617333262145e-5


class Analysis:

    def __init__(self, dataFile, temperature=300.0):
        self.dataFile = os.path.abspath(dataFile)
        self.fileDir = os.path.dirname(self.dataFile)
        self.temperature = temperature
        # self.activationEnergies = np.full(MAX_TRANSITION_STATES, 100.0, dtype=float)
        self.activationEnergies = []

    def _getActivationEnergies(self):
        with open(self.dataFile, "r") as f:
            for line in f:
                if "SPS" in line:
                    # Parse the activation energy from the line
                    energy = line.split()[-1]
                    if (float(energy) < 10.0):
                        self.activationEnergies.append(float(energy))
        f.close()

        # Get the average activation energy
        self._getAverageActivationEnergy()

    def _getAverageActivationEnergy(self):
        totalFrequency = 0
        sumEnergyFreq = 0
        for energy in self.activationEnergies:
            if energy < 10:
                freq = math.exp(-energy / (BOLTZMANN_CONSTANT * self.temperature))
                totalFrequency += freq
                sumEnergyFreq += freq * energy
        self.averageActivationEnergy = sumEnergyFreq / totalFrequency

    def _makeHistogram(self):
        mpl.style.use("default")

        histoFig, histoAxis = plt.subplots(nrows=1, ncols=1)
        histoFig.set_size_inches(10, 8)
        histoAxis.hist(
            self.activationEnergies,
            bins="auto",
            color="grey",
            alpha=0.7,
            edgecolor="black",
            linewidth=1.0,
            label="Activation Energy",
        )
        histoAxis.axvline(
            x=self.averageActivationEnergy,
            color="blue",
            linestyle="--",
            label="Average Activation Energy",
        )
        histoAxis.text(
            self.averageActivationEnergy + 0.05,
            -0.25,
            f"{self.averageActivationEnergy:.4f}",
            va="bottom",
            ha="right",
            color="blue",
        )
        histoAxis.set_xlabel("Activation Energy (eV)", fontsize=14)
        histoAxis.set_ylabel("Count", fontsize=14)
        histoAxis.set_title(
            f"Activation Energy Histogram at {self.temperature} K", fontsize=20
        )
        histoAxis.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

        histoFig.savefig(f"{self.fileDir}/histogram.png", dpi=400, transparent=False)

    def saveHistogram(self):
        # First parse the data file for activation energies
        self._getActivationEnergies()
        # Create the histogram of the activation energies
        self._makeHistogram()
