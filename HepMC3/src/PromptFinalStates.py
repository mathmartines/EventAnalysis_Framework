"""Definition of the class to identify the prompt neutrinos"""

import pyhepmc
from typing import List


class PromptFinalStates:
    """Finds prompt particles not coming from hadronic decays."""

    # PIDs of the particles that can be prompt
    _pids = [11, 12, 13, 14, 22]

    def __init__(self):
        self._memo = {}     # Cache all the particles checked in the event

    def select_prompt_particles(self, particles: List[pyhepmc.GenParticle]):
        """Excludes neutrinos from tau or hadron decays."""
        self._memo.clear()  # Clear cache before each new selection
        prompt_particles = []   # Stores only the prompt particles
        non_prompt_particles = []  # Stores all the non-prompt particles

        # Searches for prompt particles
        for particle in particles:
            if particle.abs_pid in self._pids and self.check_prompt_particle(particle):
                prompt_particles.append(particle)
            else:
                non_prompt_particles.append(particle)

        # Returns the prompt particles
        return prompt_particles, non_prompt_particles

    def check_prompt_particle(self, particle: pyhepmc.GenParticle):
        """Checks if the particle originated from a hadron decay."""
        # Particle id in the event
        part_id = particle.id

        # Check if the paricle was already checked before
        if part_id in self._memo:
            return self._memo[part_id]

        # Looks for all the incoming particles to the production vertex
        for parent in particle.production_vertex.particles_in:
            # Recursion: if the particle is not a hadron, the parents might be
            # OBS: abs(pid) > 100 -> hadrons, status == 4 -> beam particle
            if (parent.abs_pid > 100 and parent.status != 4) or not self.check_prompt_particle(parent):
                self._memo[part_id] = False
                return False

        # At this stage: the particle is prompt
        self._memo[part_id] = True
        return True

