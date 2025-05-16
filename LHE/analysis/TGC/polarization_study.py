"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import json
import copy


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WW/Longitudinal"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM"
    ]

    # Invariant mass of the diboson system
    bin_edges = [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 1000000]
    mvv_dist = ObservableHistogram(
        bin_edges=bin_edges, observable=InvariantMassObs(part_pids=[24])
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mvv_dist)

    # Book the histogram for the signal
    histograms_mvv_efts = {
        eft_term: copy.copy(mvv_dist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 10):
            # name of the .lhe file
            file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            histograms_mvv_efts[eft_term] += (xsection / number_of_evts) * current_hist

            print(f"Cross-section: {xsection} pb")
            print("Number of events per bin:")
            nevts = "\n".join([
                f"bin [{bin_edges[i]:.0f}, {bin_edges[i + 1]:.0f}]: {n:.0f}" for i, n in enumerate(current_hist)
            ])
            print(nevts)

    with open(f"{folderpath}/pp_to_WLWL.json", "w") as file_:
        simuations = {term: dist.tolist() for term, dist in histograms_mvv_efts.items()}
        json.dump(simuations, file_, indent=4)
