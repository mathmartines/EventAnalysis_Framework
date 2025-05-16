"""Definition of the total phase-space for the CMS WZ analysis."""

from EventAnalysis_Framework.HepMC3.src.DressedLeptons import LeptonsDresser
from EventAnalysis_Framework.HepMC3.src.PromptFinalStates import PromptFinalStates
import itertools


class ParticleSelectorCMS:
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


def find_ossf_pairs(leptons):
    """Finds all the OSSF pairs."""
    # All possiple pairs
    lepton_pairs = itertools.combinations(leptons, 2)
    # Selects only opposite sign and same flavor pairs
    ossf_pairs = filter(lambda pair: pair[0].pid == -pair[1].pid, lepton_pairs)

    return ossf_pairs


def assign_leptons(leptons):
    """Checks the existance of a least one opposite sign same flavor pair."""
    ossf_pairs = find_ossf_pairs(leptons)

    # Finds the pair whose mass is closest to the Z mass
    MZ_PDG = 91.1876
    sorted_pairs = sorted(
        ossf_pairs,
        key=lambda pair: abs((pair[0].momentum + pair[1].momentum).m() - MZ_PDG)
    )
    Zleptons = list(sorted_pairs[0])

    # The lepton from the W is the remaining one
    Wlepton = [particle for particle in leptons if particle not in Zleptons]

    return Zleptons, Wlepton


def total_phase_space_cuts(event):
    """Cuts that define the total phase-space region."""

    # Event must have three prompt leptons and one neutrino
    if len(event["leptons"]) != 3 and len(event["neutrinos"]) != 1:
        return False

    # Checks the invariant mass of the OSSF pair
    for ossf_pair in find_ossf_pairs(event["leptons"]):
        pair_momentum = ossf_pair[0].momentum + ossf_pair[1].momentum
        if pair_momentum.m() < 4:
            return False

    # Assing the leptons to their Z or W parent
    Zleptons, Wleptons = assign_leptons(event["leptons"])

    # Z boson momentum
    Zboson = Zleptons[0].momentum + Zleptons[1].momentum

    # Check the invariant mass cut
    return 60 < Zboson.m() < 120


def MWZ(event) -> float:
    """Computes the invariant mass of the WZ system."""
    # Charged leptons
    leptons = event["leptons"]

    # Total momentum (excluding the neutrino)
    WZ_charged_system = leptons[0].momentum + leptons[1].momentum + leptons[2].momentum

    # Neutrino momentum
    nu = event["neutrinos"][0].momentum
    nu.e = nu.pt()  # Enerygy without the longitudinal component
    nu.px = 0       # No longitudinal component

    # Invariant mass of the WZ system
    return (WZ_charged_system + nu).m()
