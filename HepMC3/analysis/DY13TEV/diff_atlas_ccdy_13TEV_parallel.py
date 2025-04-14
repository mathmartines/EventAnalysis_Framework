"""Analysis for the single differential cross-section from ATLAS (arxiv: https://arxiv.org/abs/2502.21088v1)."""

from EventAnalysis_Framework.src.Analysis import EventLoop, EventAnalysis
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.HepMC3.analysis.DY13TEV import diff_atlas_ccdy_13TEV
import pyhepmc
import numpy as np
import json
import sys


if __name__ == "__main__":
    # Folderpath filename bin index
    folderpath, eft_term, bin_index = sys.argv[1:4]

    # Event analysis that must be applied
    event_analysis = EventAnalysis(
        particles_selection=diff_atlas_ccdy_13TEV.FinalStateSelector(),
        cuts=[
            diff_atlas_ccdy_13TEV.correct_final_state,
            diff_atlas_ccdy_13TEV.lepton_cuts,
            diff_atlas_ccdy_13TEV.neutrino_cut
        ]
    )

    # Histogram for the analysis
    bin_edges = np.array([200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000])
    bin_sizes = np.array([upper - lower for lower, upper in zip(bin_edges[: -1], bin_edges[1:])])
    mt_dist = ObservableHistogram(bin_edges=bin_edges, observable=diff_atlas_ccdy_13TEV.mt_observable)

    # Performs the loop iteration
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mt_dist)

    # Reads the cross-section for the current simulated bin
    cross_section = read_xsection(f"{folderpath}/lhe_files/diff-atlas-ccdy-13TEV-{eft_term}-{bin_index}.lhe")

    # Launches the analysis on the events
    dist, evt_number = event_loop.analyse_events(
        filename=f"{folderpath}/hepmc_files/diff-atlas-ccdy-13TEV-{eft_term}-{bin_index}.hepmc",
        event_analysis=event_analysis
    )

    # Updates the distribution
    mt_dist += (cross_section / evt_number) * dist

    # Invariant mass dists
    with open(f"{eft_term}-{bin_index}.json", "w") as file_:
        diff_xsec = mt_dist / bin_sizes
        json.dump(diff_xsec.tolist(), file_, indent=4)
        print(".json file saved")
