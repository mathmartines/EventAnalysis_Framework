"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.HepMC3.analysis.TGC.ATLAS_WW_1905_04242.fiducial_phase_space import fiducial_phase_space, leading_lepton_pt, FinalStateSelector
from EventAnalysis_Framework.src.Utilities import read_xsection
import pyhepmc
import json
import sys


if __name__ == "__main__":
    folderpath, eft_term, bin_number = sys.argv[1], sys.argv[2], sys.argv[3]

    # leading lepton observable
    bin_edges_pTlead = [27, 40, 50, 60, 70, 80, 90, 100, 110, 130, 150, 175, 220, 300, 1000]
    pTlead_dist = ObservableHistogram(
        bin_edges=bin_edges_pTlead, observable=leading_lepton_pt
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[fiducial_phase_space], particles_selection=FinalStateSelector())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=pTlead_dist)

    # name of the .lhe file
    file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_number}.lhe"
    hepmc_file = f"{folderpath}/hepmc_files/{eft_term}-bin-{bin_number}.hepmc"

    # Cross-section
    xsection = read_xsection(file_name)

    # Run the analysis on the file
    hist_nevents, number_of_evts = event_loop.analyse_events(hepmc_file, event_analysis)
    hist_pt_lead = (xsection / number_of_evts) * hist_nevents

    print(hist_nevents, xsection)

    with open(f"{eft_term}-bin-{bin_number}.json", "w") as file_:
        simuations = {eft_term: hist_pt_lead.tolist()}
        json.dump(simuations, file_, indent=4)
