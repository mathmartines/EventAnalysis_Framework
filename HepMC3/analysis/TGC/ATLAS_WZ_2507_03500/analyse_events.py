
"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.HepMC3.analysis.TGC.ATLAS_WZ_2507_03500.phase_space_cuts import (ParticleSelectorATLAS,
                                                                                              fiducial_cuts,
                                                                                              transverse_mass)
from EventAnalysis_Framework.src.Utilities import read_xsection
import pyhepmc
import json
import copy

if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WZ/ATLAS_2507_03500/SMEFTsim"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM"
    ]

    # leading lepton observable
    bin_edges_mTWZ = [0, 140, 160, 180, 210, 250, 300, 400, 500, 600, 700, 900, 100000000000]
    mtWZ_hist = ObservableHistogram(
        bin_edges=bin_edges_mTWZ, observable=transverse_mass
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[fiducial_cuts], particles_selection=ParticleSelectorATLAS())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mtWZ_hist)

    # Book the histogram for the signal
    histograms_mtWZ_efts = {
        eft_term: copy.copy(mtWZ_hist) for eft_term in eft_terms
    }

    # Iterates over the terms we need to run the analysis
    for eft_term in eft_terms:
        # Iterates over the simulated bins
        for bin_index in range(1, 7):
            # name of the .lhe file
            filename = f"{eft_term}-bin-{bin_index}"
            hepmc_file = f"{folderpath}/hepmc_files/{filename}.hepmc"
            # Cross-section
            xsection = read_xsection(
                path_to_file=f"{folderpath}/banner_files/{filename}.txt",
                default_line="#  Matched Integrated weight (pb)  :"
            )
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(hepmc_file, event_analysis)
            histograms_mtWZ_efts[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

        print(histograms_mtWZ_efts[eft_term])
        print("---------------------------------------------")

    # with open(f"{folderpath}/ATLAS_WW_1905_04242_dsig_dpTlead-HepMC3.json", "w") as file_:
    #     simuations = {term: dist.tolist() for term, dist in histograms_pTlead_efts.items()}
    #     json.dump(simuations, file_, indent=4)
