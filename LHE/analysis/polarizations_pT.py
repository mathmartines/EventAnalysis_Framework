"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs, evaluate_total_momentum
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import json
import copy
import numpy as np


def pT(event: pylhe.LHEEvent) -> float:
    """Evaluates the transverse momentum of the final gauge boson."""
    w_boson = evaluate_total_momentum(event, part_pids=[24])
    return np.sqrt(np.power(w_boson[1], 2) + np.power(w_boson[2], 2))


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/PolarizationStudy/WH/AllPolarizations"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM", "CHud11", "CHud11-CHud11"
    ]

    # Transverse momentum of the final gauge bosons
    bin_edges = [75, 150, 250, 400, 600, 100000000000]
    pTV_dist = ObservableHistogram(
        bin_edges=bin_edges, observable=pT
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=pTV_dist)

    # Book the histogram for the signal
    histograms_pTV_efts = {
        eft_term: copy.copy(pTV_dist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 8):
            # name of the .lhe file
            file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            histograms_pTV_efts[eft_term] += (xsection / number_of_evts) * current_hist

            print(f"Cross-section: {xsection} pb")
            print("Number of events per bin:")
            nevts = "\n".join([
                f"bin [{bin_edges[i]:.0f}, {bin_edges[i + 1]:.0f}]: {n:.0f}" for i, n in enumerate(current_hist)
            ])
            print(nevts)

    with open(f"{folderpath}/WH-pTV.json", "w") as file_:
        simuations = {term: dist.tolist() for term, dist in histograms_pTV_efts.items()}
        json.dump(simuations, file_, indent=4)
