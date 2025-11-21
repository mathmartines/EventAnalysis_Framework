"""Collects the observables for the analysis."""

import pylhe
import numpy as np
from typing import List
from EventAnalysis_Framework.LHE.src.kinematic_funcs import evaluate_total_momentum, build_four_momentum
import abc


class FinalStateObservables(abc.ABC):
    """Abstract class to represent observables computed with a selected set of particles."""

    def __init__(self, part_pids: List[int]):
        # Absolute PIDs of the particles that must be included
        self._pids = part_pids

    @abc.abstractmethod
    def __call__(self, event: pylhe.LHEEvent):
        """Computes the observable."""
        pass


class InvariantMassObs(FinalStateObservables):
    """Computes the invariant mass for a set of particles in the event."""

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the invariant mass."""
        total_momentum = evaluate_total_momentum(event, self._pids)
        # Invariant mass
        invariant_mass_squared = np.sum(
            [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(total_momentum)]
        )
        return np.sqrt(invariant_mass_squared)


class TransverseMassObs(FinalStateObservables):
    """Computes the transverse mass of a set of particles"""

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the invariant mass."""
        # Selected particles
        particles = np.array([p for p in event.particles if abs(p.id) in self._pids])

        # Momentum of the final particles
        momentums = np.array([build_four_momentum(part) for part in particles])

        # Transverse momentum of each of the particles
        pTs = np.array([np.sqrt(p[1]**2 + p[2]**2) for p in momentums])

        # Transvese energy
        Et = np.array([np.sqrt(pt**2 + part.m**2) for pt, part in zip(pTs, particles)])

        # Transverse mass of the event
        return np.sqrt(np.power(np.sum(Et), 2) - np.power(np.sum(pTs), 2))


class TransverseMomentum(FinalStateObservables):
    """Computes the invariant mass for a set of particles in the event."""

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the invariant mass."""
        total_momentum = evaluate_total_momentum(event, self._pids)
        # Invariant mass
        transverse_mass_sq = np.power(total_momentum[1], 2) + np.power(total_momentum[2], 2)
        return np.sqrt(transverse_mass_sq)
