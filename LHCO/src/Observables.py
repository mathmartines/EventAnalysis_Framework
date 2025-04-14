"""Observalbles that can be used for the LHCO analysis"""


import vector
import numpy as np
from typing import List
from EventAnalysis_Framework.LHCO.src.EventInfo import Event


class InvariantMass:
    """Computes the invariant mass of a given set of particles"""

    def __init__(self, particles: List[str]):
        self._part_names = particles

    def _compute_inv_mass(self, event: Event):
        """Computes the invariant mass of the given event"""
        total_momentum = np.sum(vector.MomentumNumpy4D([
            particle.momentum()
            for part_name in self._part_names for particle in event.__getattr__(part_type=part_name)
        ]))
        # returns the invariant mass
        return total_momentum.m

    def __call__(self, event: Event):
        return self._compute_inv_mass(event)




