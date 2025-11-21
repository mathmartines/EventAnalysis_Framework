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
    # Folder where the files are stored
    path_to_folder = "/home/martines/work/MG5_aMC_v2_9_23/PhD/FCC-hh/pp_mumu/lhe_files"

    # Books the histogram
    bin_edges = [0, 500, 1000, 2000, 3500, 5000, 100000000000]
    invariant_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=InvariantMassObs(part_pids=[13]))

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=invariant_mass_hist)

    # Leptoquark masses
    lq_masses = [1, 4, 6, 8, 10, 20]

    # Book the histogram for the signal
    histogram_signal = {mass: copy.copy(invariant_mass_hist) for mass in lq_masses}

    # Iterates over all masses and all bins
    for mass in lq_masses:
        for bin_index in range(1, 7):
            # name of the .lhe file
            file_name = f"{path_to_folder}/mU1-{mass}TeV-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            histogram_signal[mass] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

    # Saves the json file
    # with open(f"{path_to_folder}/Mll_hist.json", "w") as file_:
    #     json.dump(histogram_signal.tolist(), file_, indent=4)
