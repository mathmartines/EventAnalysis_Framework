"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, pT, eta
import copy
import json
import numpy as np
import pylhe


def get_weight(event: pylhe.LHEEvent):
    """Get the weights of the event"""
    return event.eventinfo.weight


if __name__ == "__main__":
    # Path to the folder where the simulations are stored
    path_to_folder = "/home/martines/work/MG5_aMC_v2_9_23/PhD/HighPT-VH/ZH/lhe_files"

    # Simulated terms
    eft_terms = [
        "SM",
        # "Cphi1", "Delta4F", "CBW",  "C2JB",
        # "Cphi1-Cphi1", "Cphi1-CBW", "Cphi1-C2JB", "Cphi1-Delta4F",
        # "Delta4F-Delta4F", "Delta4F-CBW", "Delta4F-C2JB",
        # "CBW-CBW", "CBW-C2JB",
        # "C2JB-C2JB"
    ]

    # Books the histogram
    bin_edges = np.arange(300, 1600, 100)
    invariant_mass_hist = ObservableHistogram(
        bin_edges=bin_edges,
        observable=InvariantMassObs(part_pids=[25, 23]),
        get_weight=get_weight
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=invariant_mass_hist)

    # Books the histogram for each eft term
    eft_hists = {eft_term: copy.copy(invariant_mass_hist) for eft_term in eft_terms}

    # Iterates over all the terms
    for eft_term in eft_terms:
        # Iterates over all bins
        for bin_index in range(1, len(bin_edges)):
            # name of the lhe file
            file_name = f"{path_to_folder}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term (in fb)
            eft_hists[eft_term] += current_hist / number_of_evts
        print(eft_hists[eft_term])

    # Saves the json file
    with open(f"{path_to_folder}/Simulations.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)
