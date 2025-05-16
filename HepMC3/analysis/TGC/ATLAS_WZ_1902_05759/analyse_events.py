
"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.HepMC3.analysis.TGC.ATLAS_WZ_1902_05759.phase_space_cuts import ParticleSelectorATLAS, fiducial_cuts, transverse_mass
from EventAnalysis_Framework.src.Utilities import read_xsection
import pyhepmc
import json
import copy

if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WZ/CMS_2110_11231"

    # Effective terms we need to run the analysis on
    eft_terms = [
        "SM"
    ]

    # leading lepton observable
    bin_edges_mTWZ = [0, 140, 180, 250, 450, 600, 100000000000]
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
            file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_index}.lhe"
            hepmc_file = f"{folderpath}/hepmc_files/{eft_term}-bin-{bin_index}.hepmc"
            # Cross-section
            xsection = read_xsection(file_name)
            # Run the analysis on the file
            current_hist, number_of_evts = event_loop.analyse_events(hepmc_file, event_analysis)
            histograms_mtWZ_efts[eft_term] += (xsection / number_of_evts) * current_hist
            print(current_hist, xsection)

        print(histograms_mtWZ_efts[eft_term])
        print("---------------------------------------------")

    # with open(f"{folderpath}/ATLAS_WW_1905_04242_dsig_dpTlead-HepMC3.json", "w") as file_:
    #     simuations = {term: dist.tolist() for term, dist in histograms_pTlead_efts.items()}
    #     json.dump(simuations, file_, indent=4)
