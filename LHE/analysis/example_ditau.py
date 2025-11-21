"""Reads the events from a .lhe file and constructs the invariant mass distribution."""

from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.src.Utilities import read_xsection
import numpy as np
import matplotlib.pyplot as plt
import pylhe


def m(event: pylhe.LHEEvent) -> float:
    """Computes the invariant mass of the tau pair in the final state."""
    # Selects all the taus in the final state
    taus = [tau for tau in event.particles if abs(tau.id) == 15]

    # Contructs the total four-momentum of the tau pair
    # 2 x 4 matrix - each line represents the momentum of a particle, and the row their components
    taus_momentum = np.array([
        [getattr(tau, comp) for comp in "e px py pz".split()]
        for tau in taus
    ])

    # Component-wise sum -- gives the four-momentum of the tau pair
    total_momentum = np.sum(taus_momentum, axis=0)

    # Computes the invariant mass squared
    inv_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(total_momentum)]
    )

    return np.sqrt(inv_mass_squared)


if __name__ == "__main__":
    # Path where the .lhe file is stored
    lhe_file = "/Users/martines/Desktop/Physics/MG5_aMC_v2_9_20/test_tautau/Events/run_01/unweighted_events.lhe"

    # Defines distribution to construct the histogram
    bin_edges = np.arange(0, 1200, 50)
    mtautau_hist = ObservableHistogram(
        bin_edges=bin_edges,
        observable=m
    )

    # Reads the xsection
    xsec = read_xsection(path_to_file=lhe_file)
    print(f"Cross-section: {xsec} pb")

    # Total number of events
    nevents = pylhe.read_num_events(filepath=lhe_file)

    # Loops over all the events in the file
    for evt in pylhe.read_lhe(filepath=lhe_file):
        # Apply kinematical cuts
        # ...
        # Updates the histogram if passed all the cuts
        mtautau_hist.update_hist(evt)

    # Calculating the cross-section prediction for each bin
    xsec_per_bin = mtautau_hist * (xsec / nevents)

    # Plot histogram
    # Central values of the bins
    bin_central_val = np.array([np.mean([low, high]) for low, high in zip(bin_edges[:-1], bin_edges[1:])])
    plt.hist(bin_central_val, bins=bin_edges, weights=mtautau_hist, histtype='step')
    plt.xlabel("m [GeV]")
    plt.ylabel("sigma [pb]")
    plt.yscale("log")
    plt.show()
