"""Definition of the class to identify the prompt neutrinos"""

import pyhepmc
from typing import List


class PromptFinalStates:
    """Finds the prompt neutrinos or leptons that originated from a W or Z decays."""

    def select_prompt_particles(self, particles: List[pyhepmc.GenParticle]):
        """Excludes neutrinos from tau or hadron decays."""
        prompt_particles = []   # Stores only the prompt one

        # Searches for prompt particles
        for particle in particles:
            if self.check_prompt_particle(particle):
                prompt_particles.append(particle)

        # Returns the prompt particles
        return prompt_particles

    def check_prompt_particle(self, particle: pyhepmc.GenParticle):
        """Checks if the particle originated from a hadron decay."""
        # Looks for all the incoming particles to the production vertex
        for parent in particle.production_vertex.particles_in:
            # Recursion: if the particle is not a hadron, the parents might be
            if (parent.abs_pid > 100 and parent.status != 4) or not self.check_prompt_particle(parent):
                return False
        # At this stage: the particle is prompt
        return True
