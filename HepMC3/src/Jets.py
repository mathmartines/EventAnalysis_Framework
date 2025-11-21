"""Handles the construction of Jets"""

import fastjet
import pyhepmc
from typing import List


class JetsBuilder:
    """Constructs the Jets out of a list of hepmc particles."""

    def __init__(self, jet_def=fastjet.antikt_algorithm, jet_radius=0.4):
        """Parameters for the jet construction."""
        self._jet_def = fastjet.JetDefinition(jet_def, jet_radius)
        self._cluster = None    # It has to be kept in the scope if we need to access info on jets

    def cluster_particles(self, particles: List[pyhepmc.GenParticle], min_pt: float):
        """Cluster the jets"""
        # Convert HepMC3 particles to PseudoJets
        pseudo_jets = [
            fastjet.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
            for p in particles
        ]
        for jet, particle in zip(pseudo_jets, particles):
            jet.set_user_index(particle.pid)  # Store PID in user_index

        self._cluster = fastjet.ClusterSequence(pseudo_jets, self._jet_def)

        # Return jets above min_pt
        return self._cluster.inclusive_jets(min_pt)
