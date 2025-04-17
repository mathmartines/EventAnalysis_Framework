""" Implements the particle selection and cuts that are needed for the analysis."""

from EventAnalysis_Framework.LHCO.src.EventInfo import Event
from itertools import combinations
import numpy as np
import vector


def number_of_leptons(event: Event) -> bool:
    """Only three leptons in the event"""
    if len(event.electrons + event.muons) != 3:
        return False
    for lepton in event.electrons + event.muons:
        eta_cut = 2.5 if abs(lepton.typ) == 1 else 2.4
        if abs(lepton.eta) > eta_cut:
            return False
    return True


def check_ossf_pair(event: Event) -> bool:
    """Check the existance of at least one opposite sign same flavor lepton pair"""
    # Looks for a single OSSF pair
    for lepton_flavor in ["electrons", "muons"]:
        leptons_charges = [lepton.ntrk for lepton in event.__getattr__(lepton_flavor)]
        if len(leptons_charges) > 1 > abs(sum(leptons_charges)):
            return True
    return False


def missing_energy_cut(event: Event) -> bool:
    """Applies the missing energy cut."""
    if len(event.met) < 1:
        return False
    met = event.met[0]
    return met.pt > 30


def assign_leptons_decays(event: Event):
    """Finds the pair of leptons that originated from a Z boson decay and the one from the W decay"""
    # Simplest case: two leptons of one flavor and one of other
    if len(event.electrons) > 0 and len(event.muons) > 0:
        z_leptons = event.electrons if len(event.electrons) > len(event.muons) else event.muons
        w_leptons = event.muons if len(event.electrons) > len(event.muons) else event.electrons
        return z_leptons, w_leptons

    # All leptons in the event (all the same flavor at this point)
    leptons = event.electrons + event.muons

    # Stores the possible pairs and their respective invariant masses
    candidate_pairs = []

    # mass of the Z boson
    mass_z = 91.1876  # GeV

    # Finds all possible pairs
    for pair in combinations(leptons, 2):
        pair_momentum = np.sum(vector.MomentumNumpy4D([part.momentum() for part in pair]))
        candidate_pairs.append({
            "leptons": list(pair),
            "mass": abs(pair_momentum.m - mass_z)
        })

    # Finds the pair whose invariant mass is closest to the Z mass
    leptons_z = min(candidate_pairs, key=lambda cand_pair: cand_pair["mass"])["pair"]
    leptons_z = sorted(leptons_z, key=lambda particle: particle.pt, reverse=True)

    # Lepton from the W decay
    lepton_w = [lepton for lepton in leptons if lepton not in leptons_z]

    # Returns the chosen pairs
    return leptons_z, lepton_w


def leptons_cut(event: Event) -> bool:
    """Applies the cuts on the leptons."""
    # Assign the leptons
    z_leptons, w_leptons = assign_leptons_decays(event)

    # Invariant of the z pair
    inv_mass = np.sum(vector.MomentumNumpy4D([part.momentum() for part in z_leptons])).m
    if abs(inv_mass - 91.1876) > 15:
        return False

    # Leptons pT cuts
    return z_leptons[0].pt > 25 and z_leptons[1].pt > 10 and w_leptons[0].pt > 25


def btag_veto(event: Event) -> bool:
    """Veto events with a b-tagged jet"""
    return not any([jet.btag > 0 for jet in event.jets])


def min_invariant_mass(event: Event) -> bool:
    """Requirement on the min mll for any lepton pair"""
    # All leptons in the event
    leptons = event.electrons + event.muons

    # Finds all possible pairs
    for lepton_pair in combinations(leptons, 2):
        # Invariant mass of the lepton pair
        inv_mass = np.sum(vector.MomentumNumpy4D([part.momentum() for part in lepton_pair])).m
        if inv_mass < 4:
            return False

    return True


def invariant_mass_all_leptons(event: Event) -> bool:
    """Invariant mass requirement for all leptons"""
    # All leptons in the event
    leptons = event.electrons + event.muons
    leptons_momentum = np.sum(vector.MomentumNumpy4D([part.momentum() for part in leptons]))
    return leptons_momentum.m > 100
