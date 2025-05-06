"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs, evaluate_total_momentum
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import numpy as np
import copy


def invariant_mass_wz(event):
    """Computes the invariant mass of the event."""
    leptons = evaluate_total_momentum(event, [11, 13])

    # Invariant mass
    invariant_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(leptons)]
    )
    return np.sqrt(invariant_mass_squared)


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WW/CMS_2009_00119-CC/lhe_files"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "CLHud21-CLHud21"
    ]

    # Histogram for the analysis
    bin_edges = [0, 100, 200, 300, 400, 500, 600, 700, 800, 1000, 10000000000]
    mll_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=invariant_mass_wz
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mll_hist)

    # Book the histogram for the signal
    histograms_efts = {
        eft_term: copy.copy(mll_hist) for eft_term in eft_terms
    }

    total_xsec = 0

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 11):
            # name of the .lhe file
            file_name = f"{folderpath}/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            total_xsec += xsection
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            histograms_efts[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

    print(histograms_efts["SM"]/total_xsec)
