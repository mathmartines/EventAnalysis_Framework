"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.HepMC3.analysis.TGC.CMS_WZ_2110_11231.total_phase_space import ParticleSelectorCMS, total_phase_space_cuts, MWZ
from EventAnalysis_Framework.src.Utilities import read_xsection
import pyhepmc
import json
import sys


if __name__ == "__main__":
    folderpath, eft_term, bin_number = sys.argv[1], sys.argv[2], sys.argv[3]

    # leading lepton observable
    bin_edges_mWZ = [100, 160, 200, 300, 600, 3000]
    mWZ_hist = ObservableHistogram(
        bin_edges=bin_edges_mWZ, observable=MWZ
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[total_phase_space_cuts], particles_selection=ParticleSelectorCMS())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mWZ_hist)

    # name of the .lhe file
    file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_number}.lhe"
    hepmc_file = f"{folderpath}/hepmc_files/{eft_term}-bin-{bin_number}.hepmc"

    # Cross-section
    xsection = read_xsection(file_name)

    # Run the analysis on the file
    hist_nevents, number_of_evts = event_loop.analyse_events(hepmc_file, event_analysis)
    hist_mTWZ = (xsection / number_of_evts) * hist_nevents

    print(hist_nevents, xsection)

    with open(f"{eft_term}-bin-{bin_number}.json", "w") as file_:
        simuations = {eft_term: hist_mTWZ.tolist()}
        json.dump(simuations, file_, indent=4)
