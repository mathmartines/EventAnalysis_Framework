"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.HepMC3.analysis.TGC.ATLAS_WZ_1902_05759.phase_space_cuts import ParticleSelectorATLAS, fiducial_cuts, transverse_mass
from EventAnalysis_Framework.src.Utilities import read_xsection
import pyhepmc
import json
import sys


if __name__ == "__main__":
    folderpath, eft_term, bin_number = sys.argv[1], sys.argv[2], sys.argv[3]

    # leading lepton observable
    bin_edges_mTWZ = [0, 140, 180, 250, 450, 600, 100000000000]
    mtWZ_hist = ObservableHistogram(
        bin_edges=bin_edges_mTWZ, observable=transverse_mass
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[fiducial_cuts], particles_selection=ParticleSelectorATLAS())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mtWZ_hist)

    # name of the .lhe file
    file_name = f"{folderpath}/lhe_files/{eft_term}-bin-{bin_number}.lhe"
    hepmc_file = f"{folderpath}/hepmc_files/{eft_term}-bin-{bin_number}.hepmc"

    # Cross-section
    xsection = read_xsection(file_name)

    # Run the analysis on the file
    hist_nevents, number_of_evts = event_loop.analyse_events(hepmc_file, event_analysis)
    hist_mTWZ = (xsection * 1000 / number_of_evts) * hist_nevents

    print(hist_nevents, xsection)

    with open(f"{eft_term}-bin-{bin_number}.json", "w") as file_:
        simuations = {eft_term: hist_mTWZ.tolist()}
        json.dump(simuations, file_, indent=4)
