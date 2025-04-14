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
    path_to_folder = "/home/martines/work/MG5_aMC_v3_1_1/PhD/DY13TEV/PartonLevelDists/CC/lhe_files"

    eft_terms = [
        "SM", "Delta4F",
        "Delta4F-Delta4F",
        "C3psi4D2"
    ]

    # Books the histogram
    bin_edges = np.arange(0, 7200, 200)
    invariant_mass_hist = ObservableHistogram(bin_edges=bin_edges, observable=InvariantMassObs(part_pids=[11, 12]))

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=invariant_mass_hist)

    # Books the histogram for each eft term
    eft_hists = {term: copy.copy(invariant_mass_hist) for term in eft_terms}

    # Iterates over all the terms
    for eft_term in eft_terms:
        # Iterates over all bins
        for bin_index in range(1, 12):
            # name of the lhe file
            file_name = f"{path_to_folder}/parton_level_cc-{eft_term}-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            eft_hists[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

    # Saves the json file
    with open(f"{path_to_folder}/cc-parton-level.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)
