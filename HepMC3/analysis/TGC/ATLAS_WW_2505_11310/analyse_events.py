"""Event selection for ATLAS 2505.11310"""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.HepMC3.analysis.TGC.ATLAS_WW_2505_11310.fidutial_phase_space import FinalStates, event_selection
import pyhepmc
import json
import numpy as np
import sys


def mTemu(event):
    """Computes the transverse mass of the emu pair in the event."""
    prompt_leptons = event["prompt-lepton"]
    met = event["MET"][0]

    # Momentum of the electron-muon system
    leptons = prompt_leptons[0].momentum + prompt_leptons[1].momentum

    # Transverse energy of the e mu system
    et_emu = np.sqrt(np.power(leptons.m(), 2) + np.power(leptons.pt(), 2))

    return np.sqrt(np.power(met.pt() + et_emu, 2) - np.power((met + leptons).pt(), 2))


if __name__ == "__main__":
    filename = sys.argv[1]

    folderpath = '/data/01/martines/MG5_aMC_v2_9_23/PhD/TGC/WW/ATLAS_2505_11310/Validation'

    # Transverse mass of the event
    bin_edges_mTemu = [85, 200, 300, 450, 600, 1200, 1000000000]
    mTemu_hist = ObservableHistogram(bin_edges=bin_edges_mTemu, observable=mTemu)

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[event_selection], particles_selection=FinalStates())

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mTemu_hist)

    # Cross-section
    xsection = read_xsection(
        path_to_file=f"{folderpath}/banner_files/{filename}.txt",
        default_line="#  Matched Integrated weight (pb)  :"
    )

    # Performs the event analysis
    current_hist, number_of_evts = event_loop.analyse_events(
        filename=f"{folderpath}/hepmc_files/{filename}.hepmc", event_analysis=event_analysis
    )
    print(current_hist, xsection)
    print(number_of_evts)

    # Updates the histogram
    mTemu_hist = current_hist * (xsection / number_of_evts)

    with open(f"{filename}.json", "w") as file_:
        json.dump(mTemu_hist.tolist(), file_, indent=4)
