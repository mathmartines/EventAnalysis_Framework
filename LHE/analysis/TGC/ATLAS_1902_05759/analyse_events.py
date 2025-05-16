"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.LHE.analysis.TGC.ATLAS_1902_05759 import fiducial_phase_space
import json
import pylhe
import copy


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WZ/CMS_2110_11231"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM",
        # C1Hq
        "C1Hq11",
        "C1Hq12",
        "C1Hq13",
        "C1Hq22",
        "C1Hq23",
        "C1Hq33",
        "C1Hq11-C1Hq11",
        "C1Hq12-C1Hq12",
        "C1Hq13-C1Hq13",
        "C1Hq22-C1Hq22",
        "C1Hq23-C1Hq23",
        "C1Hq33-C1Hq33",
        # C3Hq
        "C3Hq11",
        "C3Hq12",
        "C3Hq13",
        "C3Hq22",
        "C3Hq23",
        "C3Hq33",
        "C3Hq11-C3Hq11",
        "C3Hq12-C3Hq12",
        "C3Hq13-C3Hq13",
        "C3Hq22-C3Hq22",
        "C3Hq23-C3Hq23",
        "C3Hq33-C3Hq33",
        # CHd
        "CHd11",
        "CHd12",
        "CHd13",
        "CHd22",
        "CHd23",
        "CHd33",
        "CHd11-CHd11",
        "CHd12-CHd12",
        "CHd13-CHd13",
        "CHd22-CHd22",
        "CHd23-CHd23",
        "CHd33-CHd33",
        # CHu
        "CHu11",
        "CHu12",
        "CHu13",
        "CHu22",
        "CHu23",
        "CHu11-CHu11",
        "CHu12-CHu12",
        "CHu13-CHu13",
        "CHu22-CHu22",
        "CHu23-CHu23",
        # CHud
        "CHud11",
        "CHud12",
        "CHud13",
        "CHud21",
        "CHud22",
        "CHud23",
        "CHud11-CHud11",
        "CHud12-CHud12",
        "CHud13-CHud13",
        "CHud21-CHud21",
        "CHud22-CHud22",
        "CHud23-CHud23"
    ]

    # Histogram for the analysis
    bin_edges = [0, 140, 180, 250, 450, 600, 100000000000]
    mtWZ_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=fiducial_phase_space.transverse_mass
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[fiducial_phase_space.fiducial_cuts])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mtWZ_hist)

    # Book the histogram for the signal
    histograms_efts = {
        eft_term: copy.copy(mtWZ_hist) for eft_term in eft_terms
    }

    total_xsec = 0

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 7):
            # name of the .lhe file
            file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            total_xsec += xsection
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            histograms_efts[eft_term] += (xsection * 1000 / number_of_evts) * current_hist
            print(current_hist, xsection)

    # Saves results
    with open(f"{folderpath}/ATLAS_WZ_1902_05759_dsig_dmTWZ.json", "w") as file_:
        simuations = {term: dist.tolist() for term, dist in histograms_efts.items()}
        json.dump(simuations, file_, indent=4)
