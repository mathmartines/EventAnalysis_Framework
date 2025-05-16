"""Fiducial phase-space definition for ATLAS WW 1905.04242."""

import pylhe
import numpy as np
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, pT, eta, M


def fiducial_phase_space(event: pylhe.LHEEvent) -> bool:
    """Cuts that define the fiducial phase space for the analysis."""
    # Leptons in the event
    leptons = [lep for lep in event.particles if abs(lep.id) in [11, 13]]
    # Neutrinos in the event
    neutrinos = [nu for nu in event.particles if abs(nu.id) in [12, 14]]
    # Missing energy
    met = np.sum([build_four_momentum(nu) for nu in neutrinos], axis=0)

    # Two final state leptons
    if len(leptons) != 2:
        return False

    # Opposite sign and opposite flavor leptons restriction
    if leptons[0].id == leptons[1].id or leptons[0].id * leptons[1].id > 0:
        return False

    # Lepton cuts
    lep1_mom = build_four_momentum(leptons[0])
    lep2_mom = build_four_momentum(leptons[1])

    if pT(lep1_mom) <= 27 or abs(eta(lep1_mom)) >= 2.5:
        return False
    if pT(lep2_mom) <= 27 or abs(eta(lep2_mom)) >= 2.5:
        return False

    # Missing energy cut
    if pT(met) <= 20:
        return False

    # Invariant mass cut of the lepton pair
    if M(lep1_mom + lep2_mom) <= 55:
        return False

    # Transverse momentum cut of the lepton pair
    if pT(lep1_mom + lep2_mom) <= 30:
        return False

    return True


def leading_lepton_pT(event: pylhe.LHEEvent) -> float:
    """Computes the momentum of the leading lepton."""
    leptons_pt = [pT(build_four_momentum(lep)) for lep in event.particles if abs(lep.id) in [11, 13]]
    return max(leptons_pt)







