"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.src.Observables import TransverseMomentum, evaluate_total_momentum
import copy
import json
import itertools
import numpy as np
import pylhe


def higgs_rapidity_cut(event: pylhe.LHEEvent):
    """Applies the cut on the rapidity of the Higgs"""
    higgs_momentum = evaluate_total_momentum(event, [25])
    # Higgs rapidity
    yhiggs = 0.5 * np.log((higgs_momentum[0] + higgs_momentum[3]) / (higgs_momentum[0] - higgs_momentum[3]))
    # Apply the cut
    return abs(yhiggs) < 2.5


if __name__ == "__main__":
    # Path to the folder where the simulations are stored
    path_to_folder = "/home/martines/work/MG5_aMC_v2_9_23/PhD/STXS/ATLAS_2410_19611/WH"

    # Simulated terms
    eft_terms = [
        "C3Hq11-CHud11"
    ]

    # Books the histogram
    bin_edges = [75, 150, 250, 400, 600, 1000000000]
    trans_momentum_v = ObservableHistogram(
        bin_edges=bin_edges,
        observable=TransverseMomentum(part_pids=[11, 12, 13, 14, 15, 16])
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[higgs_rapidity_cut])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=trans_momentum_v)

    # Books the histogram for each eft term
    eft_hists = {eft_term: copy.copy(trans_momentum_v) for eft_term in eft_terms}

    # Iterates over all the terms
    for eft_term in eft_terms:
        # Iterates over all bins
        for bin_index in range(1, 8):
            # name of the lhe file
            file_name = f"{path_to_folder}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            eft_hists[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

    # Saves the json file
    with open(f"{path_to_folder}/ATLAS_WH_2410_19611_dsigma_dpTW.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)
