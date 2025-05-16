"""Invariant mass distributions for the CC and NC channels."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
import copy
import json
import pylhe


if __name__ == "__main__":
    # Path to the folder where the simulations are stored
    path_to_folder = "/home/martines/work/MG5_aMC_v2_9_23/PhD/STXS/ATLAS_2410_19611/WH"

    # Simulated terms
    eft_terms = [
        "SM", "C3Hq11", "C3Hq11-C3Hq11", "CHud11-CHud11"
    ]

    # Books the histogram
    bin_edges = [0, 250, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 100000000000000000]
    mWh = ObservableHistogram(
        bin_edges=bin_edges,
        observable=InvariantMassObs(part_pids=[11, 12, 13, 14, 15, 16, 25])
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mWh)

    # Books the histogram for each eft term
    eft_hists = {eft_term: copy.copy(mWh) for eft_term in eft_terms}

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
            print(f"Cross-section: {xsection} pb")
            print("Number of events per bin:")
            nevts = "\n".join([
                f"bin [{bin_edges[i]:.0f}, {bin_edges[i + 1]:.0f}]: {n:.0f}" for i, n in enumerate(current_hist)
            ])
            print(nevts)

    # Saves the json file
    with open(f"{path_to_folder}/ATLAS_WH_2410_19611_dsigma_dMWH.json", "w") as file_:
        eft_hists = {term: dist.tolist() for term, dist in eft_hists.items()}
        json.dump(eft_hists, file_, indent=4)
