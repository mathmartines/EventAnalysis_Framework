"""Analysis for the single differential cross-section from ATLAS (arxiv: https://arxiv.org/abs/2502.21088v1)."""

from EventAnalysis_Framework.HepMC3.src.DressedLeptons import LeptonsDresser
from EventAnalysis_Framework.HepMC3.src.Jets import JetsBuilder
import pyhepmc


class FinalStateSelector:
    """Selects the lepton and neutrino for the analysis."""

    particles_pids = {"leptons": [11, 13], "photons": [22], "neutrinos": [12, 14, 16], "stable_final_states": []}

    def __init__(self):
        # Responsible to dress the final state leptons
        self.lepton_dresser = LeptonsDresser([11, 13], 0.1)
        self.jet_builder = JetsBuilder()

    def select_particles(self, event):
        """Selects the final state particles."""
        # Selects all the final state particles of interested
        particles = {cat: [] for cat in self.particles_pids.keys()}

        # Selects particles for the analysis
        for particle in event.particles:
            if particle.status == 1:
                particles[self.find_particle_category(particle.abs_pid)].append(particle)

        # Sort the vectors by pT - leptons and photons
        for cat in ["leptons", "photons"]:
            particles[cat] = sorted(particles[cat], key=lambda part: part.momentum.pt(), reverse=True)

        # Creates the dressed leptons - keeping only the hardest lepton
        dressed_leptons, photons = self.lepton_dresser.create_dressed_leptons(
            leptons=particles["leptons"],
            photons=particles["photons"]
        )

        # Adds the remaining photons to the list of stable final state particles
        particles["stable_final_states"].extend(photons)

        # keeps only the leptons with pT > 27 and |eta| < 2.7
        hard_leptons = []
        for dressed_lepton in dressed_leptons:
            if dressed_lepton.momentum.pt() > 27 and abs(dressed_lepton.momentum.eta()) < 2.5:
                hard_leptons.append(dressed_lepton)
            else:
                particles["stable_final_states"].append(dressed_lepton)

        # Cluster jets
        clustered_jets = self.jet_builder.cluster_particles(particles["stable_final_states"], 35)

        # Objects for the analysis
        selected_objects = {
            "leptons": hard_leptons, "neutrinos": particles["neutrinos"], "jets": clustered_jets
        }

        # Order leptons by pT again
        selected_objects["leptons"] = sorted(selected_objects["leptons"],
                                             key=lambda part: part.momentum.pt(),
                                             reverse=True)

        # All particles needed for the analysis
        return selected_objects

    def __call__(self, event):
        return self.select_particles(event)

    @classmethod
    def find_particle_category(cls, pid):
        """Finds which category the particle belongs to."""
        for category_label, pids in cls.particles_pids.items():
            if pid in pids and category_label != "stable_final_states":
                return category_label
        return "stable_final_states"


def fiducial_phase_space(event):
    """Defines the fiducial phase space for the analysis."""
    # Two leptons in the final state
    if len(event["leptons"]) != 2:
        return False

    # No jets in the event with pT >= 35 and |eta| < 4.5
    for jet in event["jets"]:
        if jet.pt() >= 35 and abs(jet.eta()) < 4.5:
            return False

    # leptons
    lead_lep = event["leptons"][0]
    sublead_lep = event["leptons"][1]
    ll_system = lead_lep.momentum + sublead_lep.momentum

    # Opposite flavor leptons
    if lead_lep.abs_pid == sublead_lep.abs_pid or lead_lep.pid * sublead_lep.pid > 0:
        return False

    # Missing energy momentum
    met = pyhepmc.FourVector(0, 0, 0, 0)
    for nu in event["neutrinos"]:
        met = met + nu.momentum

    # Missing energy cut
    if met.pt() <= 20:
        return False

    # Invariant mass cut
    if ll_system.m() <= 55:
        return False

    # Transverse momentum cut
    if ll_system.pt() <= 30:
        return False

    return True


def leading_lepton_pt(event) -> float:
    """Returns the pT of the leading lepton"""
    return event["leptons"][0].momentum.pt()
