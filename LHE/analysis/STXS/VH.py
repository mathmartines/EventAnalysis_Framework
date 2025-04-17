"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.src.Observables import TransverseMass, evaluate_total_momentum
import copy
import json
import numpy as np
import pylhe


def higgs_rapidity_cut(event: pylhe.LHEEvent):
    """Applies the cut on the rapidity of the Higgs"""
    higgs_momentum = evaluate_total_momentum(event, [25])
    # Higgs rapidity
    yhiggs = 0.5 * np.log((higgs_momentum[0] + higgs_momentum[3]) / (higgs_momentum[0] - higgs_momentum[3]))
    # Apply the cut
    return abs(yhiggs) < 2.5


def construct_eft_term_names(term: str, indices: str):
    """Constructs the name of the eft term"""
    return "-".join([f"{coef}{indices}" for coef in term.split("-")])


if __name__ == "__main__":
    path_to_folder = "/home/martines/work/MG5_aMC_v3_1_1/PhD/STXS/ZH/lhe_files"

    eft_terms = [
        "C3Hq", "C3Hq-C3Hq"
    ]
    flavor_indices = ["11", "12"]

    # Books the histogram
    # bin_edges = [75, 150, 250, 400, 600, 1000000000]
    bin_edges = [0, 250, 400, 1000, 1500, 2000, 3000, 100000000]
    trans_momentum_v = ObservableHistogram(
        bin_edges=bin_edges,
        observable=TransverseMass(part_pids=[11, 12, 13, 14, 15, 16])
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=trans_momentum_v)

    # Books the histogram for each eft term
    eft_hists = {
        construct_eft_term_names(term, index): copy.copy(trans_momentum_v)
        for term in eft_terms for index in flavor_indices
    }

    # Iterates over all the terms
    for eft_term in eft_terms:
        for flavor_index in flavor_indices:
            eft_term_name = construct_eft_term_names(eft_term, flavor_index)
            # Iterates over all bins
            for bin_index in range(1, 8):
                # name of the lhe file
                file_name = f"{path_to_folder}/{eft_term}-{flavor_index}-{bin_index}.lhe"
                # Cross-section
                xsection = read_xsection(file_name)
                # Run the analysis on the file
                current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
                # Updates the histogram for the current term
                eft_hists[eft_term_name] += (xsection / number_of_evts) * current_hist
                print(current_hist, xsection)

    # Saves the json file
    # with open(f"{path_to_folder}/STXS_ZH.json", "w") as file_:
    #     eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
    #     json.dump(eft_hists, file_, indent=4)
