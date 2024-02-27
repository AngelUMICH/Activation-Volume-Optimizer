# This script is used to plot the results of the dimer searces for
# a defect cluster or point defect. The Analysis class takes in the
# data file from OPLD and creates a histogram of the transition states.

import argparse

from analysis import Analysis

# Parse the command line for the data file name
parser = argparse.ArgumentParser()
parser.add_argument("data_file", type=str, help="data file from OPLD")
parser.add_argument("temperature", default=300.0, type=float, help="temperature in K")
args = parser.parse_args()
data_file = args.data_file
temperature = args.temperature

plotter = Analysis(data_file, temperature)
plotter.saveHistogram()
