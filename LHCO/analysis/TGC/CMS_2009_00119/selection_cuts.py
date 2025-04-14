""" Implements the particle selection and cuts that are needed for the analysis."""


from EventAnalysis_Framework.LHCO.src.EventInfo import Event
import numpy as np
import vector


def select_objects(event: Event):
    """Selects all the objects needed for the event analysis."""
    # Particles for the analysis
    selected_particles = Event([])

    # Electrons for the event
    for electron in event.electrons:
        if electron.pt > 10 and (abs(electron.eta) < 1.479 or 1.566 < abs(electron.eta) < 2.5):
            selected_particles.append(electron)

    # Muons for the event
    for muon in event.muons:
        if muon.pt > 10 and abs(muon.eta) < 2.4:
            selected_particles.append(muon)

    # Missing energy
    selected_particles.extend(event.met)

    # Jets
    for jet in event.jets:
        # jets tagged as b-jets must be included - later we must veto events with b-tagged jets
        if jet.btag > 0:
            selected_particles.append(jet)
        else:
            # pT of the jets must be greater than 30 GeV and |eta| < 4.7
            if jet.pt < 30 or abs(jet.eta) > 4.7:
                continue

            # jet momentum
            jet_momentum = jet.momentum()

            # Only include jets with DeltaR(j, l) > 0.4
            include_jet = True
            for lepton in selected_particles.electrons + selected_particles.muons:
                if jet_momentum.deltaR(lepton.momentum()) < 0.4:
                    include_jet = False
                    break

            if include_jet:
                selected_particles.append(jet)

    return selected_particles


def opposite_sign_lepton_pair(event: Event) -> bool:
    """Selects only events with two opposite sign and opposite flavor leptons."""
    # Only two leptons in the final state
    if len(event.electrons + event.muons) != 2:
        return False

    # Opposite flavor leptons only
    if len(event.electrons) != len(event.muons):
        return False

    # Opposite sign leptons
    return event.electrons[0].ntrk * event.muons[0].ntrk < 0


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




