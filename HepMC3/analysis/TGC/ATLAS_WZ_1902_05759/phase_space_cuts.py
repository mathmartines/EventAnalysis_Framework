"""Definition of the phase-space accorfing to ATLAS-1902.05759."""

from EventAnalysis_Framework.HepMC3.src.DressedLeptons import LeptonsDresser
from EventAnalysis_Framework.HepMC3.src.PromptFinalStates import PromptFinalStates
import pyhepmc
import numpy as np

# Masses and widths from PDG
MZ_PDG = 91.1876
MW_PDG = 80.385
GammaZ_PDG = 2.4952
GammaW_PDG = 2.085


class ParticleSelectorATLAS:
    """Selects the lepton and neutrino for the analysis."""

    particles_pids = {"leptons": [11, 13], "photons": [22], "neutrinos": [12, 14]}
    pids_for_analysis = [11, 13, 22, 12, 14]

    def __init__(self):
        # Responsible to dress the final state leptons
        self.lepton_dresser = LeptonsDresser([11, 13], 0.1)
        self.prompt_final_states = PromptFinalStates()

    def select_particles(self, event):
        """Selects the final state particles."""
        # Selects all the final state particles of interested
        particles = {cat: [] for cat in self.particles_pids.keys()}

        # Selects particles needed for the analysis
        for particle in event.particles:
            part_pid = particle.abs_pid
            if particle.status == 1 and part_pid in self.pids_for_analysis:
                part_cat = self.find_particle_category(part_pid)
                particles[part_cat].append(particle)

        # Creates the dressed leptons
        dressed_leptons, photons = self.lepton_dresser.create_dressed_leptons(
            leptons=particles["leptons"],
            photons=particles["photons"]
        )

        # Finds the prompt leptons and prompt neutrinos
        prompt_leptons = self.prompt_final_states.select_prompt_particles(dressed_leptons)
        prompt_nu = self.prompt_final_states.select_prompt_particles(particles["neutrinos"])

        # Objects for the analysis
        selected_objects = {
            "leptons": prompt_leptons, "neutrinos": prompt_nu
        }

        # Order by pT
        selected_objects = {
            cat: sorted(parts, key=lambda part: part.momentum.pt(), reverse=True)
            for cat, parts in selected_objects.items()
        }

        # All particles needed for the analysis
        return selected_objects

    def __call__(self, event):
        return self.select_particles(event)

    @classmethod
    def find_particle_category(cls, pid):
        """Finds which category the particle belongs to."""
        for category_label, pids in cls.particles_pids.items():
            if pid in pids:
                return category_label


# !!!!!!!!!!!! Cuts !!!!!!!!!!!!!!!!

def propagator(mll, mass, width):
    """Computes the propagator square"""
    return 1. / ((mll * mll - mass * mass)**2 + (mass * width)**2)


def resonant_shape_algorithm(leptons, neutrinos):
    """Assigns the Z leptons and the W lepton-neutrino pair."""
    candidates = None
    estimator = 0

    # Neutrinos in the event
    nuW = neutrinos[0]

    # Creates all possible pairs of OSSF leptons
    for index_l1, lepZ1 in enumerate(leptons[:-1]):
        for lepZ2 in leptons[index_l1 + 1:]:
            # Check if they are OSSF
            if lepZ1.pid != -lepZ2.pid:
                continue

            # Remaining leptons
            remaining = [lep for lep in leptons if lep not in [lepZ1, lepZ2]]
            if not remaining:
                continue
            lepW = remaining[0]

            # Check if it's a valid W candidate
            if lepW.abs_pid + 1 != nuW.abs_pid or lepW.pid * nuW.pid > 0:
                continue

            # Build momentas
            Zboson = lepZ1.momentum + lepZ2.momentum
            Wboson = lepW.momentum + nuW.momentum

            # Compute estimator
            Zprop = propagator(Zboson.m(), MZ_PDG, GammaZ_PDG)
            Wprop = propagator(Wboson.m(), MW_PDG, GammaW_PDG)
            curr_estimator = Zprop * Wprop

            if curr_estimator > estimator:
                candidates = ([lepZ1, lepZ2], [lepW, nuW])
                estimator = curr_estimator

    return candidates


def fiducial_cuts(event) -> bool:
    """Returns true if the event is inside the fiducial phase space definition."""

    # All leptons and neutrinos in the event
    leptons = event["leptons"]
    neutrinos = event["neutrinos"]

    if len(leptons) != 3 or len(neutrinos) != 1:
        return False

    candidates = resonant_shape_algorithm(leptons, neutrinos)
    if candidates is None:
        return False

    z_leptons, w_leptons = candidates

    # Build 4-vectors
    Zlep1 = z_leptons[0].momentum
    Zlep2 = z_leptons[1].momentum
    Wlep = w_leptons[0].momentum
    Wnu = w_leptons[1].momentum
    Zboson = Zlep1 + Zlep2

    # W transverse mass
    norm = Wlep.pt() * Wnu.pt()
    if norm == 0:
        return False
    cos_lnu = np.cos(pyhepmc.delta_phi(Wlep, Wnu))
    Wboson_mT = np.sqrt(2 * norm * (1 - cos_lnu))

    # ----------------- Cuts -----------------
    if Zlep1.pt() <= 15 or Zlep2.pt() <= 15 or Wlep.pt() <= 20:
        return False
    if Zlep1.abs_eta() >= 2.5 or Zlep2.abs_eta() >= 2.5 or Wlep.abs_eta() >= 2.5:
        return False
    if abs(Zboson.m() - MZ_PDG) > 10:
        return False
    if Wboson_mT <= 30:
        return False
    if pyhepmc.delta_r_eta(Zlep1, Zlep2) <= 0.2:
        return False
    if pyhepmc.delta_r_eta(Zlep1, Wlep) <= 0.3 or pyhepmc.delta_r_eta(Zlep2, Wlep) <= 0.3:
        return False

    return True


def transverse_mass(event) -> float:
    """Computes the transverse mass of the event."""
    particles = event["leptons"] + event["neutrinos"]
    # Total pT
    WZ_pT = np.sum([particle.momentum.pt() for particle in particles])
    # px
    WZ_px = np.sum([particle.momentum.px for particle in particles])
    WZ_py = np.sum([particle.momentum.py for particle in particles])
    # Transverse mass of the event
    return np.sqrt(WZ_pT**2 - WZ_px**2 - WZ_py**2)
