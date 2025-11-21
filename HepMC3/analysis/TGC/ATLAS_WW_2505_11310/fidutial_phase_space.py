"""Defines the fiducial phase-space region where the measurement is performed."""
import fastjet
import pyhepmc
from EventAnalysis_Framework.HepMC3.src.DressedLeptons import LeptonsDresser
from EventAnalysis_Framework.HepMC3.src.PromptFinalStates import PromptFinalStates
from EventAnalysis_Framework.HepMC3.src.Jets import JetsBuilder
from collections import defaultdict
from typing import Optional, List, Callable


def classify_particles(particles, classifier: Callable):
    """Classify the particles according to classifier."""
    classified_particles = defaultdict(list)
    for particle in particles:
        cat = classifier(particle)
        if cat:
            classified_particles[cat].append(particle)
    return classified_particles


class FinalStates:
    """Selects the final state particles and objects for the analysis."""

    # Absolute PIDs for the final state particles we need
    _particle_pids = {
        "lepton": [11, 13],
        "neutrino": [12, 14],
        "photon": [22],
    }

    # Cuts for the prompt leptons
    _prompt_pt_cut = 27
    _loose_pt_cut = 10
    _abs_eta_cut = 2.5

    def __init__(self):
        """Initializes all the objects needed to find the final state particles."""
        self._prompt_part_selector = PromptFinalStates()
        self._leptons_dresser = LeptonsDresser(delta_r=0.1)
        self._jet_builder = JetsBuilder()

    def select_particles(self, event: pyhepmc.GenEvent):
        """Selects all the particles that are needed for the event analysis."""
        # Select only final state particles
        final_particles = [p for p in event.particles if p.status == 1]

        # Select prompt particles (excluding hadron/tau decays)
        prompt_particles, non_prompt_particles = self._prompt_part_selector.select_prompt_particles(final_particles)

        # Organize prompt particles into categories
        particles = classify_particles(prompt_particles, self.find_category)

        # Dress the leptons
        dressed_leptons, remaining_photons = self._leptons_dresser.create_dressed_leptons(
            leptons=particles["lepton"],
            photons=particles["photon"]
        )

        # Classify leptons (prompt-lepton or loose-lepton)
        classified_leptons = classify_particles(dressed_leptons, self.classify_lepton)
        for cat, leps in classified_leptons.items():
            particles[cat].extend(leps)

        # Use remaining photons for jets
        jet_inputs = non_prompt_particles + remaining_photons

        # Cluster jets
        jets = self._jet_builder.cluster_particles(jet_inputs, min_pt=20)
        for jet in jets:
            jet.set_user_index(0)   # No tagged jets at this point

        # Identify b-quarks from hard process (for b-tagging) - Only one hard scatter quark
        bquark = [p for p in event.particles if p.abs_pid == 5 and p.status == 23]
        # Apply the tagging
        if len(bquark):
            self.btag_jets(jets, bquark[0])

        # Classify jets (jet or b-jet)
        classified_jets = classify_particles(jets, self.classify_jet)
        for cat, jets in classified_jets.items():
            particles[cat].extend(jets)

        # Removes the keys that are not needed anymore
        particles.pop("lepton")
        particles.pop("photon")

        # Constructing the missing energy
        met_vis = pyhepmc.FourVector(0, 0, 0, 0)

        # From leptons
        for cat in ["prompt-lepton", "loose-lepton"]:
            for particle in particles[cat]:
                met_vis = met_vis - particle.momentum

        # From jets
        for cat in ["jet", "b-jet"]:
            for jet in particles[cat]:
                jet_mom = pyhepmc.FourVector(jet.px(), jet.py(), jet.pz(), jet.E())
                met_vis = met_vis - jet_mom

        # Adds the MET to the list of particles
        particles["MET"] = [met_vis]

        return dict(particles)

    @staticmethod
    def btag_jets(jets: List[fastjet.PseudoJet], parton: pyhepmc.GenParticle, dr_threshold: float = 0.4) -> List[fastjet.PseudoJet]:
        """
        Tags jets that are within ΔR < dr_threshold of a b-parton.

        Only the closest jet (if within threshold) is tagged with user_index = 1.
        Others are marked as 0 (assumed that all jets are already tagged with 0).
        """
        parton_mom = parton.momentum
        closest_jet = None
        min_dr = float('inf')

        # Looks for the closest jet to the parton
        for jet in jets:
            jet_mom = pyhepmc.FourVector(jet.px(), jet.py(), jet.pz(), jet.E())
            dr = pyhepmc.delta_r_eta(parton_mom, jet_mom)
            if dr < min_dr:
                min_dr = dr
                closest_jet = jet

        # Apply tagging
        if min_dr < dr_threshold:
            closest_jet.set_user_index(1)

        return jets

    @classmethod
    def find_category(cls, particle: pyhepmc.GenParticle) -> Optional[str]:
        """Find the particle category given its absolute PID."""
        for category, pid_list in cls._particle_pids.items():
            if particle.abs_pid in pid_list:
                return category
        return None

    @classmethod
    def classify_lepton(cls, lepton: pyhepmc.GenParticle) -> Optional[str]:
        """Classify a dressed lepton as prompt or loose."""
        pt = lepton.momentum.pt()
        eta = lepton.momentum.abs_eta()

        if eta > cls._abs_eta_cut:
            return None
        elif pt > cls._prompt_pt_cut:
            return "prompt-lepton"
        elif pt > cls._loose_pt_cut:
            return "loose-lepton"
        return None

    @classmethod
    def classify_jet(cls, jet: fastjet.PseudoJet) -> Optional[str]:
        """Classify the jet as jet of b-jet."""
        pt = jet.pt()
        eta = abs(jet.eta())
        btag = jet.user_index()

        if btag and pt > 20 and eta < 2.5:
            return "b-jet"
        elif pt > 30 and eta < 4.5:
            return "jet"
        return None

    def __call__(self, event: pyhepmc.GenEvent):
        return self.select_particles(event)


def event_selection(event) -> bool:
    """Applies all the cuts. Returns True if the event passed all cuts, False otherwise."""

    prompt_leptons = event["prompt-lepton"]
    loose_leptons = event["loose-lepton"]
    b_jets = event["b-jet"]

    # Two opposite sign different flavor leptons
    if len(prompt_leptons) != 2:
        return False
    if prompt_leptons[0].abs_pid == prompt_leptons[1].abs_pid or prompt_leptons[0].pid * prompt_leptons[1].pid > 0:
        return False

    # No additional loose-lepton in the event
    if len(loose_leptons):
        return False

    # No b-jets
    if len(b_jets):
        return False

    # Invariant mass of the electron-muon system
    leptons_mom = prompt_leptons[0].momentum + prompt_leptons[1].momentum

    return leptons_mom.m() > 85

