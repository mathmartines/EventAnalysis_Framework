"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.LHE.src.Observables import evaluate_total_momentum
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHE.analysis.TGC.CMS_2110_11231 import analysis_functions
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import numpy as np
import copy


def invariant_mass_wz(event):
    """Computes the invariant mass of the event."""
    leptons = evaluate_total_momentum(event, [11, 13])

    neutrinos = np.array([
        [getattr(part, comp) for comp in "e px py pz".split()]
        for part in event.particles if abs(part.id) in [12, 14]
    ])

    for index, momentum in enumerate(neutrinos):
        transv_momentum = np.sqrt(np.power(momentum[1], 2) + np.power(momentum[2], 2))
        neutrinos[index][0] = transv_momentum
        neutrinos[index][-1] = 0

    total_neutrinos = np.sum(neutrinos, axis=0)

    total_momentum = leptons + total_neutrinos

    # Invariant mass
    invariant_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(total_momentum)]
    )
    return np.sqrt(invariant_mass_squared)


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WZ/CMS_2110_11231_no_jets/lhe_files"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM"
    ]

    # Histogram for the analysis
    bin_edges = [100, 160, 200, 300, 600, 3000]
    mll_hist = ObservableHistogram(
        bin_edges=bin_edges, observable=invariant_mass_wz
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[
        analysis_functions.number_of_leptons,
        analysis_functions.min_invariant_mass,
        analysis_functions.invariant_mass_Zleptons

    ])  # no cuts being applied

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mll_hist)

    # Book the histogram for the signal
    histograms_efts = {
        eft_term: copy.copy(mll_hist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 8):
            # name of the .lhe file
            file_name = f"{folderpath}/{eft_term}-bin-{bin_index}.lhe"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(file_name, event_analysis)
            # Updates the histogram for the current term
            histograms_efts[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

        total_xsec = np.sum(histograms_efts[eft_term])

        print(histograms_efts[eft_term], total_xsec)
        print(histograms_efts[eft_term]/total_xsec)
