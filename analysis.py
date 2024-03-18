import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

BOLTZMANN_CONSTANT = 8.617333262145e-5


class Analysis:

    def __init__(self, dataDirectory, temperature=300.0):
        self.dataDirectory = os.path.abspath(dataDirectory)
        self.temperature = temperature
        self.activationEnergies = []
        self.averageActivationEnergies = []
        self.atomicVolumeRadii = []

    def _getAvRadius(self, filename):
        return float(os.path.splitext(filename)[0].split("_")[-1])

    def _getActivationEnergies(self, filename):
        self.atomicVolumeRadii.append(self._getAvRadius(filename))
        energies = []
        with open(filename, "r") as f:
            for line in f:
                if "SPS" in line:
                    # Parse the activation energy from the line
                    energy = line.split()[-1]
                    if float(energy) < 10.0:
                        energies.append(float(energy))
        f.close()
        self.activationEnergies.append(energies)
        e = self._getAverageActivationEnergy(energies)
        self.averageActivationEnergies.append(e)

    def _getAverageActivationEnergy(self, energies):
        totalFrequency = 0
        sumEnergyFreq = 0
        for energy in energies:
            if energy < 10:
                freq = math.exp(-energy / (BOLTZMANN_CONSTANT * self.temperature))
                totalFrequency += freq
                sumEnergyFreq += freq * energy
        return sumEnergyFreq / totalFrequency

    def _makeHistogram(self):
        mpl.style.use("default")
        clr = ["navy", "orange", "grey"]
        dataLabels = [f"AV radius {self.atomicVolumeRadii[i]:.2f} | Avg Energy {average:.4f} eV" for i,average in enumerate(self.averageActivationEnergies)]

        histoFig, histoAxis = plt.subplots(nrows=1, ncols=1)
        histoFig.set_size_inches(10, 8)
        histoAxis.hist(
            self.activationEnergies,
            bins="auto",
            color=clr[:3],
            alpha=0.7,
            edgecolor="black",
            linewidth=2.0,
            label=dataLabels
        )
        for index, average in enumerate(self.averageActivationEnergies):
            histoAxis.axvline(
                x=average,
                color=clr[index],
                linestyle="--",
            )
        histoAxis.set_xlabel("Activation Energy (eV)", fontsize=14)
        histoAxis.set_ylabel("Count", fontsize=14)
        histoAxis.set_title(
            f"Activation Energy Histogram at {self.temperature} K", fontsize=20
        )
        histoAxis.legend()
        histoAxis.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

        histoFig.savefig(f"{self.dataDirectory}/histogram.png", dpi=200, transparent=False)

    def saveHistogram(self):
        # First parse the data file for activation energies for each file
        for file in os.listdir(self.dataDirectory):
            if file.endswith(".log"):
                self._getActivationEnergies(self.dataDirectory+file)
        # Create the histogram of the activation energies
        self._makeHistogram()
