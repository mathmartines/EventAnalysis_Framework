""" Implements the particle selection and cuts that are needed for the analysis."""

from EventAnalysis_Framework.LHCO.src.EventInfo import Event
import numpy as np
import vector


def number_of_leptons(event: Event) -> bool:
    """Only three leptons in the event"""
    return len(event.electrons + event.muons) == 3


def check_ossf_pair(event: Event) -> bool:
    """Check the existance of at least one OSSF pair"""



def leptons_pt_cuts(event: Event) -> bool:
    """Applies the pT cuts on the selected leptons."""
    electron_pt = event.electrons[0].pt
    muon_pt = event.muons[0].pt

    # Leading and subleading pTs
    lead_pt = electron_pt if electron_pt > muon_pt else muon_pt
    sublead_pt = muon_pt if electron_pt > muon_pt else electron_pt

    # Leading and subleating cuts
    return lead_pt > 25 and sublead_pt > 20


def missing_energy_cut(event: Event) -> bool:
    """Applies the missing energy cut."""
    if len(event.met) < 1:
        return False
    met = event.met[0]
    return met.pt > 20


def lepton_pair_cuts(event: Event) -> bool:
    """Applies the cut on mll and pTll"""
    electron = event.electrons[0]
    muon = event.muons[0]
    # Total momentum of the lepton pair
    leptons_momentum = electron.momentum() + muon.momentum()
    # Cuts
    return leptons_momentum.m > 20 and leptons_momentum.pt > 30


def btag_veto(event: Event) -> bool:
    """Veto events with a b-tagged jet"""
    return not any([jet.btag > 0 for jet in event.jets])


def number_of_jets(event: Event) -> bool:
    """Veto events with more than 1 jet"""
    return len(event.jets) < 2


def projected_ptmiss_cut(event: Event) -> bool:
    """Cut on the projected pTmiss"""
    # Missing et vector
    met_momentum = event.met[0].momentum()

    # Find the closest lepton to the met momentum
    leptons = event.electrons + event.muons

    # Closest momentum
    closest_mom = leptons[0].momentum()
    deltaR = met_momentum.deltaR(closest_mom)

    for lepton in leptons[1:]:
        lepton_momentum = lepton.momentum()
        # Current distance
        curr_dist = met_momentum.deltaR(lepton_momentum)
        # Distance from the met
        if curr_dist < deltaR:
            deltaR = curr_dist
            closest_mom = lepton_momentum

    # Computes the azimuthal angle difference
    if abs(met_momentum.deltaphi(closest_mom)) > np.pi/2:
        return True

    # Calculates the orthogonal component of pT miss with respect to pT
    met_3vector = met_momentum.to_3D()

    # trasverse momentum of the lepton
    lepton_pt = vector.MomentumNumpy3D(
        (closest_mom.px, closest_mom.py, 0), dtype=[("px", float), ("py", float), ("pz", float)])
    lepton_pt_unit = lepton_pt.unit()

    # Perpendicular component
    met_3vector_perp = met_3vector - met_3vector.dot(lepton_pt_unit) * lepton_pt_unit

    # Applies the cut
    return met_3vector_perp.mag > 20




