"""Cuts that define the fiducial phase space region."""

import pylhe
import numpy as np
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum, pT, eta, M, deltaR

# Masses and widths from PDG
MZ_PDG = 91.1876
MW_PDG = 80.385
GammaZ_PDG = 2.4952
GammaW_PDG = 2.085


def propagator(mll, mass, width):
    """Computes the propagator square"""
    return 1. / ((mll * mll - mass * mass)**2 + (mass * width)**2)


def resonant_shape_algorithm(leptons, neutrinos):
    """Assigns the Z leptons and the W lepton-neutrino pair."""
    candidates = None
    estimator = 0

    # Neutrino in the event
    nuW = neutrinos[0]

    # Creates all possible pairs of OSSF leptons
    for index_l1, lepZ1 in enumerate(leptons[:-1]):
        for lepZ2 in leptons[index_l1 + 1:]:
            # Check if they are OSSF
            if lepZ1.id != -lepZ2.id:
                continue

            # Remaining lepton
            remaining = [lep for lep in leptons if lep not in [lepZ1, lepZ2]]
            if not remaining:
                continue
            lepW = remaining[0]

            # Check if it's a valid W candidate
            if abs(lepW.id) + 1 != abs(nuW.id) or lepW.id * nuW.id > 0:
                continue

            # Build momenta
            Zboson = build_four_momentum(lepZ1) + build_four_momentum(lepZ2)
            Wboson = build_four_momentum(lepW) + build_four_momentum(nuW)

            # Compute estimator
            Zprop = propagator(M(Zboson), MZ_PDG, GammaZ_PDG)
            Wprop = propagator(M(Wboson), MW_PDG, GammaW_PDG)
            curr_estimator = Zprop * Wprop

            if curr_estimator > estimator:
                candidates = ([lepZ1, lepZ2], [lepW, nuW])
                estimator = curr_estimator

    return candidates


def fiducial_cuts(event: pylhe.LHEEvent) -> bool:
    """Returns true if the event is inside the fiducial phase space definition."""

    # All leptons and neutrinos in the event
    leptons = [lep for lep in event.particles if abs(lep.id) in [11, 13]]
    neutrinos = [nu for nu in event.particles if abs(nu.id) in [12, 14]]

    if len(leptons) != 3 or len(neutrinos) != 1:
        return False

    candidates = resonant_shape_algorithm(leptons, neutrinos)
    if candidates is None:
        return False

    z_leptons, w_leptons = candidates

    # Build 4-vectors
    Zlep1 = build_four_momentum(z_leptons[0])
    Zlep2 = build_four_momentum(z_leptons[1])
    Wlep = build_four_momentum(w_leptons[0])
    Wnu = build_four_momentum(w_leptons[1])
    Zboson = Zlep1 + Zlep2

    # Transverse mass of W
    norm = pT(Wlep) * pT(Wnu)
    if norm == 0:
        return False
    cos_lnu = np.dot(Wlep[1:3], Wnu[1:3]) / norm
    Wboson_mT = np.sqrt(2 * norm * (1 - cos_lnu))

    # ----------------- Cuts -----------------
    if pT(Zlep1) <= 15 or pT(Zlep2) <= 15 or pT(Wlep) <= 20:
        return False
    if abs(eta(Zlep1)) >= 2.5 or abs(eta(Zlep2)) >= 2.5 or abs(eta(Wlep)) >= 2.5:
        return False
    if abs(M(Zboson) - MZ_PDG) > 10:
        return False
    if Wboson_mT <= 30:
        return False
    if deltaR(Zlep1, Zlep2) <= 0.2:
        return False
    if deltaR(Zlep1, Wlep) <= 0.3 or deltaR(Zlep2, Wlep) <= 0.3:
        return False

    return True


def transverse_mass(event: pylhe.LHEEvent) -> float:
    """Computes the transverse mass of the event."""
    particles_4vec = [build_four_momentum(part) for part in event.particles if abs(part.id) in [11, 12, 13, 14]]

    WZ_pT = np.sum([pT(vec) for vec in particles_4vec])
    WZ_px, WZ_py = np.sum(particles_4vec, axis=0)[1:3]

    return np.sqrt(WZ_pT**2 - WZ_px**2 - WZ_py**2)
