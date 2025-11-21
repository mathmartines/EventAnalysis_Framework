"""
    Reads all the events from the file and builds invariant mass histograms.
    The .lhe file contains events with different weights.
"""
import json
import pylhe
import numpy as np
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, eta, pT
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop


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
    # Path to where the files are store
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/USMEFT-HLLHC/DY/Zprime/lhe_files_10_TeV"

    # Books the histogram
    bin_edges = [500, 750, 1000, 1250, 1500, 2000, 3000, 4000, 1000000]
    inv_mass_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=InvariantMassObs(part_pids=[11]),
        get_weight=get_weight
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[leptons_kin_cuts])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=inv_mass_hist)

    # Iterates over all bins
    for bin_index in range(1, 6):
        # name of the lhe file
        file_name = f"{folderpath}/Zprime-beta-0E-01-MXX-10-TeV-bin-{bin_index}.lhe"
        # Run the analysis on the file
        current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
        # Updates the histogram for the current term (in fb)
        inv_mass_hist += 1e3 * current_hist / number_of_evts
        print(current_hist)

    # Saves the json file
    with open(f"{folderpath}/SM.json", "w") as file_:
        sm_pred = {"SM": inv_mass_hist.tolist()}
        json.dump(sm_pred, file_, indent=4)
