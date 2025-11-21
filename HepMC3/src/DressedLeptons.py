"""
    Creates the dressed leptons by adding the four-momenta of the photons within a distance R
    from the leptons.
"""

from typing import List
import pyhepmc


def select_particles(event, pids):
    """Returns the list of sorted particles."""
    return sorted([part for part in event if part.abs_pid in pids],
                  key=lambda particle: particle.momentum.pt(), reverse=True)


class LeptonsDresser:
    """Creates the dressed leptons for the event analysis."""

    def __init__(self, delta_r: float):
        """
        :param delta_r: size of the radius around the leptons that any photon momenta should be added
                        to the lepton momenta.
        """
        self._delta_r = delta_r

    def create_dressed_leptons(self, leptons, photons):
        """Dress all the leptons in the event."""
        # List with all dressed leptons
        dressed_leptons = []

        # Dress all leptons in the event
        for lepton in leptons:
            dressed_lepton, photons = self.dress_lepton(lepton, photons)
            dressed_leptons.append(dressed_lepton)

        # Returns a list with all the dressed leptons in the event
        return dressed_leptons, photons

    def dress_lepton(self, lepton, photons):
        """
        Dress the leptons with all the photons around it.
        Returns the dressed lepton and the remaining list of photons.
        """
        # Store the photons that have been used
        used_photons = []

        # lepton direction
        dressed_momentum = lepton.momentum

        # Iterates over all photons
        for photon in photons:
            # Checks the distance
            if pyhepmc.delta_r_eta(lepton.momentum, photon.momentum) < self._delta_r:
                # Updates the lepton momentum and saves the photon for exclusion
                dressed_momentum = dressed_momentum + photon.momentum
                used_photons.append(photon)

        # Excludes the used photons
        photons = [photon for photon in photons if photon not in used_photons]
        lepton.momentum = dressed_momentum

        return lepton, photons

