"""Analysis functions for the channels CMS ...."""

import itertools
import numpy as np


# Cuts:
# Number of leptons - ok
# Identify the Z and W leptons - ok
# 60 < M(lz1, lz2) < 120 - ok
# min mll (always same flavor) < 4 GeV - ok

def invariant_mass(particles):
    """Calculates the invariant mass of a set of particles."""
    momentums = np.array([
        [getattr(part, comp) for comp in "e px py pz".split()] for part in particles
    ])
    total_momentum = np.sum(momentums, axis=0)
    # Invariant mass
    invariant_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(total_momentum)]
    )
    return np.sqrt(invariant_mass_squared)


def pt(particle):
    return np.sqrt(particle.px**2 + particle.py**2)


def number_of_leptons(event) -> bool:
    """Checks the number of leptons in the event."""
    leptons = [particle for particle in event.particles if abs(particle.id) in [11, 13, 15]]
    return len(leptons) == 3


def find_ossf_pairs(event):
    """Finds the ossf pairs"""
    leptons_pids = [11, 13, 15]

    # Search for one pair
    leptons = [particle for particle in event.particles if abs(particle.id) in leptons_pids]

    # All possiple pairs
    lepton_pairs = itertools.combinations(leptons, 2)
    # Selects only opposite sign and same flavor pairs
    ossf_pairs = filter(lambda pair: pair[0].id == -pair[1].id, lepton_pairs)

    return ossf_pairs


def assign_leptons(event):
    """Checks the existance of a least one opposite sign same flavor pair."""
    ossf_pairs = find_ossf_pairs(event)

    # Search for one pair
    leptons = [particle for particle in event.particles if abs(particle.id) in [11, 13, 15]]

    # Finds the pair whose mass is closest to the Z mass
    mZ = 9.118760e+01
    sorted_pairs = sorted(ossf_pairs, key=lambda pair: abs(invariant_mass(pair) - mZ))
    leptonsZ = list(sorted_pairs[0])

    # The lepton from the W is the remaining one
    leptonW = [particle for particle in leptons if particle not in leptonsZ]

    return leptonsZ, leptonW


def invariant_mass_Zleptons(event) -> bool:
    """Computes the invariant mass of the Z lepton pairs."""
    z_leptons, w_leptons = assign_leptons(event)
    # z_leptons = sorted(z_leptons, key=lambda part: np.sqrt(part.px**2 + part.py**2), reverse=True)

    # if pt(z_leptons[0]) < 25 or pt(z_leptons[1]) < 10 or pt(w_leptons[0]) < 25:
    #     return False

    # invariant mass check
    return 60. < invariant_mass(z_leptons) < 120.


def min_invariant_mass(event) -> bool:
    """Minimum invariant mass of any OSSF pair."""
    for pair in find_ossf_pairs(event):
        if invariant_mass(pair) < 4.:
            return False
    return True


def z_boson_pt(event):
    """Calculates the transverse momentum of the Z boson pT."""
    z_leptons, _ = assign_leptons(event)
    momentums = np.array([
        [getattr(part, comp) for comp in "e px py pz".split()] for part in z_leptons
    ])
    total_momentum = np.sum(momentums, axis=0)
    return np.sqrt(np.power(total_momentum[1], 2) + np.power(total_momentum[2], 2))