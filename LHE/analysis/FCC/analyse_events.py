"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
import copy
import json
import numpy as np
import pylhe

if __name__ == "__main__":
    path_to_folder = "/home/martines/work/MG5_aMC_v3_1_1/PhD/FCC-hh/Leptoquark-tautau/betaL/lhe_files/mU1_4TeV"

    # Books the histogram
    bin_edges = list(range(0, 16400, 400)) + [100000000000000]
    invariant_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=InvariantMassObs(part_pids=[11, 13]))

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=invariant_mass_hist)

    # Book the histogram for the signal
    histogram_signal = copy.copy(invariant_mass_hist)

    # Iterates over all bins
    for bin_index in range(1, 42):
        # name of the .lhe file
        file_name = f"{path_to_folder}/x1L-x1L-x1L-x1L-bin-{bin_index}.lhe"
        # Cross-section
        xsection = read_xsection(file_name)
        # Run the analysis on the file
        current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
        # Updates the histogram for the current term
        histogram_signal += (xsection / number_of_evts) * current_hist
        print(current_hist, xsection)

    # Saves the json file
    with open(f"{path_to_folder}/Mll_hist.json", "w") as file_:
        json.dump(histogram_signal.tolist(), file_, indent=4)
