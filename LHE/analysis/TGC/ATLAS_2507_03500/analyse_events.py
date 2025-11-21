"""Constructs the invariant mass distribution for the WW production"""

from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, pT
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Analysis import EventAnalysis, EventLoop
import pylhe
import numpy as np


def mTWZ(event: pylhe.LHEEvent):
    """Transverse momentum of the WZ system."""
    # Leptons + Neutrinos in the event
    wz_particles = np.array([build_four_momentum(part) for part in event.particles if abs(part.id) in [11, 12, 13, 14]])

    # Transverse momentum sum
    ptsum = np.sum([pT(momentum) for momentum in wz_particles])
    pxsum = np.sum([momentum[1] for momentum in wz_particles])
    pysum = np.sum([momentum[2] for momentum in wz_particles])

    # Transverse mass
    return np.sqrt(ptsum**2 - pxsum**2 - pysum**2)


if __name__ == "__main__":
    # Path to the folder where the .lhe files are stored
    folderpath = "/home/martines/work/MG5_aMC_v2_9_23/PhD/TGC/WZ/ATLAS_2507_03500/SMEFTsim/SM/Events/run_01"

    # leading lepton observable
    bin_edges_mTemu = [0, 500, 1000, 1000000000]
    mTWZ_hist = ObservableHistogram(
        bin_edges=bin_edges_mTemu, observable=mTWZ
    )

    # Constructs the event analysis
    event_analysis = EventAnalysis(cuts=[])

    # Performs the loop over the events
    event_loop = EventLoop(file_reader=pylhe.read_lhe, histogram=mTWZ_hist)

    # Performs the event analysis
    current_hist, number_of_evts = event_loop.analyse_events(
        filename=f"{folderpath}/unweighted_events.lhe", event_analysis=event_analysis
    )

    print(current_hist)



