"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.LHE.src.Observables import InvariantMassObs
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import copy


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v3_1_1/PhD/TGC/WZ/CMS_2110_11231/lhe_files"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM"
    ]

    # Histogram for the analysis
    bin_edges = [0, 100, 200, 300, 400, 500, 1000, 1500, 3000, 10000000]
    mll_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=InvariantMassObs([11, 12, 13, 14])
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mll_hist)

    # Book the histogram for the signal
    histograms_efts = {
        eft_term: copy.copy(mll_hist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 10):
            # name of the .lhe file
            file_name = f"{folderpath}/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            histograms_efts[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

    print(histograms_efts)
