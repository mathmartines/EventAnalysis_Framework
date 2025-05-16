"""Handles the construction of Jets"""

import fastjet
import pyhepmc
from typing import List


class JetsBuilder:
    """Constructs the Jets out of a list of hepmc particles."""

    def __init__(self, jet_def=fastjet.antikt_algorithm, jet_radius=0.4):
        """Parameters for the jet construction."""
        self._jet_def = fastjet.JetDefinition(jet_def, jet_radius)

    def cluster_particles(self, particles: List[pyhepmc.GenParticle], min_pt: float):
        """Cluster the jets"""
        pseudo_jets = [
            fastjet.PseudoJet(part.momentum.px, part.momentum.py, part.momentum.pz, part.momentum.e)
            for part in particles
        ]
        cluster = fastjet.ClusterSequence(pseudo_jets, self._jet_def)
        # returns the final jets
        return cluster.inclusive_jets(min_pt)
