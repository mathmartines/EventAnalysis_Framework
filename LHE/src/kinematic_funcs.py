"""Kinematic functions"""

import pylhe
import numpy as np


def build_four_momentum(particle: pylhe.LHEParticle):
    """Returns the four-momentum"""
    return np.array([getattr(particle, comp) for comp in "e px py pz".split()])


def pT(momentum) -> float:
    """Returns the pT of the particle."""
    return np.sqrt(momentum[1]**2 + momentum[2]**2)


def eta(momentum) -> float:
    """Computes the pseudo-rapidity of the particle."""
    # Norm of the 3-momentum
    p = np.sqrt(np.sum([comp ** 2 for comp in momentum[1:]]))
    # pseudo-rapidity
    return -1 / 2 * np.log((p - momentum[3]) / (p + momentum[3]))


def phi(momentum) -> float:
    """Azimuthal angle from -pi to pi"""
    px = momentum[1]
    py = momentum[2]
    return np.arctan2(py, px)


def delta_phi(p1, p2) -> float:
    """Returns the azimuthal angle difference Δφ between two momenta, in range [-π, π]."""
    dphi = phi(p1) - phi(p2)
    if dphi > np.pi:
        return dphi - 2 * np.pi
    if dphi < -np.pi:
        return dphi + 2 * np.pi
    return dphi


def delta_eta(p1, p2) -> float:
    """Returns the pseudo-rapidity difference"""
    return eta(p1) - eta(p2)


def deltaR(p1, p2) -> float:
    """Computes the DeltaR"""
    return np.sqrt(delta_eta(p1, p2)**2 + delta_phi(p1, p2)**2)


def M(momentum):
    """Computes the invariant mass."""
    # Invariant mass
    invariant_mass_squared = np.sum(
        [(1 if index == 0 else -1) * np.power(value, 2) for index, value in enumerate(momentum)]
    )
    return np.sqrt(invariant_mass_squared)


