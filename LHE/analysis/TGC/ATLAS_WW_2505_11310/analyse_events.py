"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, M, pT
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
import pylhe
import copy
import numpy as np


def mTemu(event: pylhe.LHEEvent):
    """Computes the transverse mass of the emu pair in the event."""
    # Neutrinos in the event
    neutrinos = np.sum([build_four_momentum(nu) for nu in event.particles if abs(nu.id) in [12, 14]], axis=0)
    # Charged leptons
    leptons = np.sum([build_four_momentum(lepton) for lepton in event.particles if abs(lepton.id) in [11, 13]], axis=0)

    # Transverse energy of the e mu system
    et_emu = np.sqrt(np.power(M(leptons), 2) + np.power(pT(leptons), 2))

    return np.sqrt(np.power(pT(neutrinos) + et_emu, 2) - np.power(pT(neutrinos + leptons), 2))


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WW/ATLAS_2505_11310/Validation/lhe_files"

    # EFT terms we need to check
    eft_terms = ["SM"]

    # leading lepton observable
    bin_edges_mTemu = [0, 85, 200, 300, 450, 600, 1200,  1000000000]
    mTemu_hist = ObservableHistogram(
        bin_edges=bin_edges_mTemu, observable=mTemu
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mTemu_hist)

    # Constructs the distribution for each of the EFT terms
    eft_terms_hists = {term: copy.copy(mTemu_hist) for term in eft_terms}

    # Loops over all the terms we need to launch the analysis on
    for eft_term in eft_terms:
        for bin_index in range(1, 8):
            lhe_filename = f"{eft_term}-{bin_index}.lhe"
            # cross-section
            xsec = read_xsection(path_to_file=f"{folderpath}/{lhe_filename}")
            # Performs the event analysis
            current_hist, number_of_evts = event_loop.analyse_events(
                filename=f"{folderpath}/{lhe_filename}", event_analysis=event_analysis
            )
            print(f"bin - {bin_index}:")
            print(current_hist)
            # Update hist
            eft_terms_hists[eft_term] += current_hist * (xsec / number_of_evts)

    print(eft_terms_hists)



