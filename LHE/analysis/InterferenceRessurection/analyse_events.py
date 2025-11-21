"""Constructs the asymmetry to resurect the interference terms"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.LHE.analysis.InterferenceRessurection.SpecialFrame import SpecialFrame
import pylhe
import numpy as np
import json

if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/InterferenceRessurection/gamma_jj/lhe_files"

    # EFT terms simulated
    eft_terms = ["SM"]

    # Launches the analysis over all events
    event_analysis = EventAnalysis(cuts=[])

    # Histogram for the analysis
    bin_edges = [0, np.pi/4, np.pi/2, 3 * np.pi/4, np.pi]
    phi_hist = ObservableHistogram(bin_edges=bin_edges, observable=SpecialFrame())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=phi_hist)

    # Stores the histogram for each EFT term
    hists = {}

    # Launches the analysis over the events
    for eft_term in eft_terms:
        # name of the .lhe file
        file_name = f"{folderpath}/{eft_term}.lhe"
        # Run the analysis
        hists[eft_term], number_of_evts = event_loop.analyse_events(file_name, event_analysis)
        hists[eft_term] /= number_of_evts

    with open(f"{folderpath}/gamma_jj_v2.json", "w") as file_:
        simuations = {term: dist.tolist() for term, dist in hists.items()}
        json.dump(simuations, file_, indent=4)

