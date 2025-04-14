"""Analysis for the single differential cross-section from ATLAS (arxiv: https://arxiv.org/abs/2502.21088v1)."""
import copy
from EventAnalysis_Framework.src.Analysis import EventLoop, EventAnalysis
from EventAnalysis_Framework.src.Utilities import read_xsection
from EventAnalysis_Framework.src.Histogram import ObservableHistogram
from EventAnalysis_Framework.HepMC3.src.DressedLeptons import LeptonsDresser
import pyhepmc
import numpy as np
import json


class FinalStateSelector:
    """Selects the lepton and neutrino for the analysis."""

    particles_pids = {"Leptons": [11, 13], "Photons": [22], "Neutrinos": [12, 14]}
    pids_for_analysis = [11, 13, 22, 12, 14]

    def __init__(self):
        # Responsible to dress the final state leptons
        self.lepton_dresser = LeptonsDresser([11, 13], 0.1)

    def select_particles(self, event):
        """Selects the final state particles."""
        # Selects all the final state particles of interested
        particles = {cat: [] for cat in self.particles_pids.keys()}

        for particle in event.particles:
            if particle.status == 1 and particle.abs_pid in self.pids_for_analysis:
                particles[self.find_particle_category(particle.abs_pid)].append(particle)

        # Sort the vectors by pT
        particles = {
            cat: sorted(parts, key=lambda part: part.momentum.pt(), reverse=True)
            for cat, parts in particles.items()
        }

        # Creates the dressed leptons - keeping only the hardest lepton
        dressed_leptons = self.lepton_dresser.create_dressed_leptons(
            leptons=particles["Leptons"], photons=particles["Photons"]
        )[:1]
        neutrinos = particles["Neutrinos"][:1]

        # All particles needed for the analysis
        return dressed_leptons + neutrinos

    @classmethod
    def find_particle_category(cls, pid):
        """Finds which category the particle belongs to."""
        for category_label, pids in cls.particles_pids.items():
            if pid in pids:
                return category_label


def correct_final_state(event):
    """Checks if we have the lepton and its corrected neutrino"""
    if len(event) < 2:
        return False
    lepton, neutrino = event[:2]
    # Checks the pids
    if abs(lepton.abs_pid - neutrino.abs_pid) > 1:
        return False

    return lepton.pid * neutrino.pid < 0


def lepton_cuts(event):
    """Performs the cuts on the leptons"""
    lepton = event[0]
    return lepton.momentum.pt() > 65 and abs(lepton.momentum.eta()) < 2.4


def neutrino_cut(event):
    """Performs the cuts on the neutrinos"""
    neutrino = event[1]
    return neutrino.momentum.pt() > 85


def mt_observable(event):
    """Computes the transverse mass of the event."""
    lepton, neutrino = event[:2]
    dphi = pyhepmc.delta_phi(lepton.momentum, neutrino.momentum)
    lepton_pt = lepton.momentum.pt()
    neutrino_pt = neutrino.momentum.pt()
    # transverse mass of the event
    return np.sqrt(2 * lepton_pt * neutrino_pt * (1 - np.cos(dphi)))


if __name__ == "__main__":
    # Event analysis that must be applied
    event_analysis = EventAnalysis(
        particles_selection=FinalStateSelector(),
        cuts=[correct_final_state, lepton_cuts, neutrino_cut]
    )

    # Histogram for the analysis
    bin_edges = np.array([200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000])
    bin_sizes = np.array([upper - lower for lower, upper in zip(bin_edges[: -1], bin_edges[1:])])
    mt_dist = ObservableHistogram(bin_edges=bin_edges, observable=mt_observable)

    # Performs the loop iteration
    event_loop = EventLoop(file_reader=pyhepmc.open, histogram=mt_dist)

    # folder where the simulations are stored
    folderpath = '/home/martines/work/MG5_aMC_v3_1_1/PhD/DY13TEV/diff-atlas-ccdy-13TEV/Validation/wminus'

    # EFT terms we need to launch the analysis
    eft_terms = [
        "SM", "cHl3", "cHl3-cHl3",
        "cll1", "cll1-cll1",
        "clq3", "clq3-clq3"
    ]

    # dict to store the predictions for each eft term
    eft_predictions = {term: copy.copy(mt_dist) for term in eft_terms}

    # Iterates over all the terms and bins
    for eft_term in eft_terms:
        for bin_index in range(1, 11):
            # Reads the cross-section for the current simulated bin
            cross_section = read_xsection(f"{folderpath}/lhe_files/diff-atlas-ccdy-13TEV-{eft_term}-{bin_index}.lhe")

            # Launches the analysis on the events
            dist, evt_number = event_loop.analyse_events(
                filename=f"{folderpath}/hepmc_files/diff-atlas-ccdy-13TEV-{eft_term}-{bin_index}.hepmc",
                event_analysis=event_analysis
            )
            print(dist)
            print(cross_section)

            eft_predictions[eft_term] += (cross_section / evt_number) * dist

    # Invariant mass dists
    with open('wminus.json', "w") as file_:
        simulations_hists = {term: (dist / bin_sizes).tolist() for term, dist in eft_predictions.items()}
        json.dump(simulations_hists, file_, indent=4)
        print(".json file saved")
