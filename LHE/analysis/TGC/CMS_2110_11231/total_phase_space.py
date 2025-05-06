"""Analysis functions for the channels CMS ...."""

import itertools
import numpy as np
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, M, pT
import pylhe


# Cuts:
# Number of leptons - ok
# Identify the Z and W leptons - ok
# 60 < M(lz1, lz2) < 120 - ok
# min mll (always same flavor) < 4 GeV - ok

def find_ossf_pairs(leptons):
    """Finds all the OSSF pairs."""
    # All possiple pairs
    lepton_pairs = itertools.combinations(leptons, 2)
    # Selects only opposite sign and same flavor pairs
    ossf_pairs = filter(lambda pair: pair[0].id == -pair[1].id, lepton_pairs)

    return ossf_pairs


def assign_leptons(leptons):
    """Checks the existance of a least one opposite sign same flavor pair."""
    ossf_pairs = find_ossf_pairs(leptons)

    # Finds the pair whose mass is closest to the Z mass
    MZ_PDG = 91.1876
    sorted_pairs = sorted(
        ossf_pairs,
        key=lambda pair: abs(M(build_four_momentum(pair[0]) + build_four_momentum(pair[1])) - MZ_PDG)
    )
    Zleptons = list(sorted_pairs[0])

    # The lepton from the W is the remaining one
    Wlepton = [particle for particle in leptons if particle not in Zleptons]

    return Zleptons, Wlepton


def total_phase_space_cuts(event: pylhe.LHEEvent) -> bool:
    """Applies the cuts that defines the total phase space region."""
    # Leptons in the event
    leptons = [particle for particle in event.particles if abs(particle.id) in [11, 13]]
    # Neutrinos
    neutrinos = [particle for particle in event.particles if abs(particle.id) in [12, 14]]

    # Check the number of particles in the final state
    if len(leptons) != 3 or len(neutrinos) != 1:
        return False

    # Checks the invariant mass of the OSSF pair
    for ossf_pair in find_ossf_pairs(leptons):
        pair_momentum = build_four_momentum(ossf_pair[0]) + build_four_momentum(ossf_pair[1])
        if M(pair_momentum) < 4:
            return False

    # Assing the leptons to their Z or W parent
    Zleptons, Wleptons = assign_leptons(leptons)

    # Z boson momentum
    Zboson = build_four_momentum(Zleptons[0]) + build_four_momentum(Zleptons[1])

    # Check the invariant mass cut
    return 60 < M(Zboson) < 120


def pTZ(event: pylhe.LHEEvent) -> float:
    """Computes the transverse momentum of the Z boson pair."""
    leptons = [particle for particle in event.particles if abs(particle.id) in [11, 13]]
    Zleptons, _ = assign_leptons(leptons)
    # Z boson momentum
    Zboson = build_four_momentum(Zleptons[0]) + build_four_momentum(Zleptons[1])
    # Computes the transverse momentum
    return pT(Zboson)


def MWZ(event: pylhe.LHEEvent) -> float:
    """Computes the invariant mass of the WZ system."""
    leptons = [particle for particle in event.particles if abs(particle.id) in [11, 13]]
    neutrinos = [particle for particle in event.particles if abs(particle.id) in [12, 14]]

    # Total momentum (excluding the neutrino)
    WZ_charged_system = np.sum([build_four_momentum(lepton) for lepton in leptons], axis=0)
    # Neutrino momentum
    nu = build_four_momentum(neutrinos[0])
    nu[0] = pT(nu)  # Enerygy without the longitudinal component
    nu[3] = 0       # No longitudinal component

    # Invariant mass of the WZ system
    return M(WZ_charged_system + nu)
