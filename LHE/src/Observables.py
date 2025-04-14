"""Collects the observables for the analysis."""

import pylhe
import numpy as np
from typing import List


def evaluate_total_momentum(event: pylhe.LHEEvent, part_pids: list):
    """Calculates the total momentum taking into account only the particles with PIDs in the part_pid list."""
    # Four-momentum
    momentum = np.array([
        [getattr(part, comp) for comp in "e px py pz".split()]
        for part in event.particles if abs(part.id) in part_pids
    ])
    # Total four-momentum of the event
    return np.sum(momentum, axis=0)


class InvariantMassObs:
    """Computes the invariant mass for a set of particles in the event."""

    def __init__(self, part_pids: List[int]):
        # Absolute PIDs of the particles that must be included
        self._pids = part_pids

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the invariant mass."""
        total_momentum = evaluate_total_momentum(event, self._pids)
        # Invariant mass
        invariant_mass_squared = np.sum(
            [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(total_momentum)]
        )
        return np.sqrt(invariant_mass_squared)


class TransverseMomentum:
    """Computes the invariant mass for a set of particles in the event."""

    def __init__(self, part_pids: List[int]):
        # Absolute PIDs of the particles that must be included
        self._pids = part_pids

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the invariant mass."""
        total_momentum = evaluate_total_momentum(event, self._pids)
        # Invariant mass
        transverse_mass_sq = np.power(total_momentum[1], 2) + np.power(total_momentum[2], 2)
        return np.sqrt(transverse_mass_sq)
