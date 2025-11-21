"""
    Reads all the events from the file and builds transverse mass histograms.
    The .lhe file contains events with different weights.
"""
import json
import pylhe
import numpy as np
from EventAnalysis_Framework.LHE.src.read_lhe import read_lhe
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, eta, pT
from EventAnalysis_Framework.LHE.src.Observables import TransverseMassObs
from EventAnalysis_Framework.src.Histogram import WeightedHistogramManager
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop


def leptons_kin_cuts(event: pylhe.LHEEvent):
    """Applies the cuts on the transverse momentum of the final state leptons."""
    # Momentum of the final state leptons
    leptons_4v = np.array([build_four_momentum(lepton) for lepton in event.particles if abs(lepton.id) == 11])
    # Cuts on the pT and pseudo-rapidity
    leptons_kin = all([pT(vector) > 20 and abs(eta(vector)) < 2.5 for vector in leptons_4v])
    # Veto event if needed
    return leptons_kin


def get_weights(event: pylhe.LHEEvent):
    """Get the weights of the event"""
    return event.weights


def build_filename(parameters: tuple):
    """Constructs the names of the files we need to store the pseudo-data."""
    coup, mass = parameters
    return f"CCDY-Zprime-beta_{coup * 10:.0f}E-01-MXX_{mass}-TeV.json"


if __name__ == "__main__":
    # Path to where the files are store
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/USMEFT-HLLHC/CCDY/Zprime"

    # Folders that contains the events
    folders = [f"lhe_files_{MXX}_TeV" for MXX in [2, 3, 4, 5, 6]]

    # Respective value of the coupling that was simulated in each of
    mass_coup_dict = {2: 0.3, 3: 0.5, 4: 0.5, 5: 0.8, 6: 0.8}

    # Dictionary with the values of the couplings for each reweight
    reweight_labels = {
        # Mass: labels
        # 2 TeV
        2:
            {"rwgt_1": {"beta": 0.3, "MXX": 2}, "beta_-3E-1": {"beta": -0.3, "MXX": 2}},
        3:
            {
                "rwgt_1": {"beta": 0.5, "MXX": 3}, "beta_3E-1": {"beta": 0.3, "MXX": 3},
                "beta_-3E-1": {"beta": -0.3, "MXX": 3}, "beta_-5E-1": {"beta": -0.5, "MXX": 3}
            },
        4:
            {
                "rwgt_1": {"beta": 0.5, "MXX": 4}, "beta_3E-1": {"beta": 0.3, "MXX": 4},
                "beta_8E-1": {"beta": 0.8, "MXX": 4}, "beta_-3E-1": {"beta": -0.3, "MXX": 4},
                "beta_-5E-1": {"beta": -0.5, "MXX": 4}, "beta_-8E-1": {"beta": -0.8, "MXX": 4}
            },
        5:
            {
                "rwgt_1": {"beta": 0.8, "MXX": 5}, "beta_3E-1": {"beta": 0.3, "MXX": 5},
                "beta_5E-1": {"beta": 0.5, "MXX": 5}, "beta_10E-1": {"beta": 1.0, "MXX": 5},
                "beta_-3E-1": {"beta": -0.3, "MXX": 5}, "beta_-5E-1": {"beta": -0.5, "MXX": 5},
                "beta_-83E-1": {"beta": -0.8, "MXX": 5}, "beta_-10E-1": {"beta": -1.0, "MXX": 5}
            },
        6:
            {
                "rwgt_1": {"beta": 0.8, "MXX": 6}, "beta_3E-1": {"beta": 0.3, "MXX": 6},
                "beta_5E-1": {"beta": 0.5, "MXX": 6}, "beta_10E-1": {"beta": 1.0, "MXX": 6},
                "beta_-3E-1": {"beta": -0.3, "MXX": 6}, "beta_-5E-1": {"beta": -0.5, "MXX": 6},
                "beta_-83E-1": {"beta": -0.8, "MXX": 6}, "beta_-10E-1": {"beta": -1.0, "MXX": 6},
                "beta_12E-1": {"beta": 1.2, "MXX": 6}, "beta_-12E-1": {"beta": -1.2, "MXX": 6}
            }
    }

    # Number of bins in each of the files
    nbins = 5

    # Bin edges for the histogram
    bin_edges = [500, 750, 1000, 1250, 1500, 2000, 3000, 4000, 1000000]

    # Constructs the event analysis with the required kinematic cuts
    event_analysis = EventAnalysis(cuts=[leptons_kin_cuts])

    # Initializes the histogram for each benchmark
    histograms_model = {}

    # Loops over the folders
    for foldername, (MXX, beta) in zip(folders, mass_coup_dict.items()):
        # Initializes an empty histogram for each parameter point
        for model_param in reweight_labels[MXX].values():
            histograms_model[(model_param["beta"], model_param["MXX"])] = np.zeros(len(bin_edges) - 1)

        # Prefix of the lhe file
        lhe_prefix = f"Zprime-beta-{beta * 10:.0f}E-01-MXX-{MXX}-TeV-bin"

        # Book the histogram
        histograms_mxx = WeightedHistogramManager(
            bin_edges=bin_edges, observale=TransverseMassObs(part_pids=[11, 12]), get_weights=get_weights,
            hist_names=list(reweight_labels[MXX].keys())
        )

        # Performs the loop over the events
        event_loop = EventLoop(file_reader=read_lhe, histogram=histograms_mxx)

        # Launches the event analysis in each bin of the simulation
        for bin_number in range(1, nbins + 1):
            lhe_filename = f"{folderpath}/{foldername}/{lhe_prefix}-{bin_number}.lhe"

            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(lhe_filename, event_analysis)

            # Updates the histograms
            for hist_name, model_param in reweight_labels[MXX].items():
                param_values = (model_param["beta"], model_param["MXX"])
                histograms_model[param_values] += 1e3 * current_hist[hist_name] / number_of_evts

    # Saves each of the histogram
    for model_param, hist in histograms_model.items():
        pseudo_data = {"data": hist.tolist()}
        with open(f"{folderpath}/json_files/{build_filename(model_param)}", "w") as file_:
            json.dump(pseudo_data, file_, indent=4)
