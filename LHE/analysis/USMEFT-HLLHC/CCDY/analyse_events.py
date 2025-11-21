"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.src.Observables import TransverseMassObs
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, pT, eta
import copy
import json
import numpy as np
import pylhe


def leptons_kin_cuts(event: pylhe.LHEEvent):
    """Applies the cuts on the transverse momentum of the final state leptons."""
    # Momentum of the final state leptons
    leptons_4v = np.array([build_four_momentum(lepton) for lepton in event.particles if abs(lepton.id) == 11])
    # Cuts on the pT and pseudo-rapidity
    leptons_kin = all([pT(vector) > 20 and abs(eta(vector)) < 2.5 for vector in leptons_4v])
    # Veto event if needed
    return leptons_kin


def get_weight(event: pylhe.LHEEvent):
    """Get the weights of the event"""
    return event.eventinfo.weight


if __name__ == "__main__":
    # Path to the folder where the simulations are stored
    path_to_folder = "/home/martines/work/MG5_aMC_v2_9_23/PhD/USMEFT-HLLHC/CCDY"

    # Simulated terms
    eft_terms = [
        "SM", "Cphi1", "Delta4F", "CBW",
        "Cphi1-Cphi1", "Cphi1-CBW", "Cphi1-Delta4F",
        "Delta4F-Delta4F", "Delta4F-CBW", "CBW-CBW"
    ]

    # Books the histogram
    bin_edges = [500, 750, 1000, 1250, 1500, 2000, 3000, 3500, 1000000]
    transverse_mass_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=TransverseMassObs(part_pids=[11, 12]),
        get_weight=get_weight
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[leptons_kin_cuts])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=transverse_mass_hist)

    # Books the histogram for each eft term
    eft_hists = {eft_term: copy.copy(transverse_mass_hist) for eft_term in eft_terms}

    # Iterates over all the terms
    for eft_term in eft_terms:
        # Iterates over all bins
        for bin_index in range(1, 9):
            # name of the lhe file
            file_name = f"{path_to_folder}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term (in fb)
            eft_hists[eft_term] += 1e3 * current_hist / number_of_evts
            print(current_hist, xsection)
        print(eft_hists[eft_term])
        
    # Saves the json file
    with open(f"{path_to_folder}/CCDY.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)

    with open(f"{path_to_folder}/CCDY_pseudo_data_Zprime.json", "w") as file_:
        eft_hists = {"data": eft_hists["SM"], "SM": eft_hists["SM"], "bin-edges": bin_edges}
        json.dump(eft_hists, file_, indent=4)
