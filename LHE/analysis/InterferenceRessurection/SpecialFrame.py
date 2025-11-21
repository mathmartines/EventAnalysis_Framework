"""Reconstructs the special reference frame to measure the azimuthal decay angle."""

import pylhe
import numpy as np
import vector
from EventAnalysis_Framework.LHE.src.kinematic_funcs import build_four_momentum
import abc


class LorentzTransformation(abc.ABC):
    """Abstract class to represent Lorentz transformations."""
    @abc.abstractmethod
    def apply_transformation(self, four_vector: np.array):
        """Applies the Lorentz transformation to the 4-vector."""
        pass


class LorentzBoost(LorentzTransformation):
    """Applies a Lorentz boost along the z-direction."""

    def __init__(self, beta):
        self._beta = beta
        self._gamma = np.power(1 - self._beta**2, -0.5)

    def apply_transformation(self, four_vector: np.array):
        """Applies the boost along the z-direction."""
        # E and pz in the new frame
        energy = self._gamma * (four_vector[0] - self._beta * four_vector[-1])
        pz = self._gamma * (-self._beta * four_vector[0] + four_vector[-1])
        # Returns the new four-vector
        return np.array([energy, *four_vector[1:-1], pz], dtype=np.float32)


class LorentzRotationZ(LorentzTransformation):
    """Rotation around the z-axis"""

    def __init__(self, angle: float):
        self._angle = angle

    def apply_transformation(self, four_vector: np.array):
        """Rotation around the z-axis by an angle theta."""
        px = np.cos(self._angle) * four_vector[1] - np.sin(self._angle) * four_vector[2]
        py = np.sin(self._angle) * four_vector[1] + np.cos(self._angle) * four_vector[2]
        return np.array([four_vector[0], px, py, four_vector[-1]], dtype=np.float32)


class LorentzRotationY(LorentzTransformation):
    """Rotation around the z-axis"""

    def __init__(self, angle: float):
        self._angle = angle

    def apply_transformation(self, four_vector: np.array):
        """Rotation around the z-axis by an angle theta."""
        px = np.cos(self._angle) * four_vector[1] + np.sin(self._angle) * four_vector[3]
        pz = -np.sin(self._angle) * four_vector[1] + np.cos(self._angle) * four_vector[3]
        return np.array([four_vector[0], px, four_vector[2], pz], dtype=np.float32)


class CompositeLorentzTransformations(LorentzTransformation):
    """Applies a sequence of lorentz transformations."""

    def __init__(self, transf):
        self._lorentz = transf

    def apply_transformation(self, four_vector: np.array):
        """Applies all lorentz transformations"""
        for lorentz_transf in self._lorentz:
            four_vector = lorentz_transf.apply_transformation(four_vector)
        return four_vector


class Particles:
    """Stores the information of the final state particles in the event."""
    # PIDs of the particles in each category
    _pids = {
        "gamma": [22],
        "jet": list(range(1, 6)) + [21]
    }

    def __init__(self, event: pylhe.LHEEvent):
        # Stores the final state particles of the event
        self._final_particles = {}
        # Stores the 4-vectors of all final state particles
        self._four_vectors = {}
        # Intializes all containers
        self._select_final_states(event)
        self._build_four_vectors()

    def _select_final_states(self, event: pylhe.LHEEvent):
        """Saves the final state particles."""
        # Reset all final particles
        self._final_particles = {typ: [] for typ in self._pids}

        # Search over the final state particles in the event and allocate it in its category
        for particle in event.particles:
            # Final state particles only
            if particle.status == 1:
                for typ, pids in self._pids.items():
                    if abs(particle.id) in pids:
                        self._final_particles[typ].append(particle)

    def _build_four_vectors(self):
        """Creates the four-vector of all particles involved."""
        for typ, particles in self._final_particles.items():
            self._four_vectors[typ] = np.array([build_four_momentum(p) for p in particles], dtype=np.float32)

    def get_four_vectors(self, typ: str):
        """Returns the four vectors of the particles from the category of typ."""
        return self._four_vectors.get(typ, np.array([]))

    def get_particles(self, typ: str):
        """Returns the particles from the category of typ."""
        return self._final_particles.get(typ, np.array([]))

    def apply_lorentz_transf(self, lorentz_tranf: LorentzTransformation):
        """Applies the lorentz transformations to all the particles."""
        for typ in self._final_particles:
            self._four_vectors[typ] = np.array([
                lorentz_tranf.apply_transformation(four_vec) for four_vec in self._four_vectors[typ]
            ])


class SpecialFrame:
    """Computes the four-vectors in the special reference frame."""

    def __call__(self, event: pylhe.LHEEvent):
        """Computes the azimuthal angle."""
        # Particles with four-momentum in the special reference frame
        particles = self.find_special_frame(event)

        # Finds the anti-particle
        jet1, jet2 = particles.get_particles("jet")

        # Antiparticle index
        index_right = 0 if jet1.id < 0 else 1

        # momentum of the right-handed particle
        momentum = particles.get_four_vectors("jet")[index_right]

        event_weight = event.eventinfo.weight
        # Azimuthal angle between -pi to pi and the weight of the event
        return np.abs(np.arctan2(momentum[2], momentum[1])), event_weight

    @staticmethod
    def find_special_frame(event: pylhe.LHEEvent):
        """Particles momentum in the special reference frame."""
        # Selects all the final state particles
        particles = Particles(event)

        # 4-vector of the Diboson system
        wboson_4v = np.sum(particles.get_four_vectors("jet"), axis=0)
        gamma_4v = np.sum(particles.get_four_vectors("gamma"), axis=0)
        db_4vector = wboson_4v + gamma_4v

        # Boost to go to the CM frame of the DB system
        beta_z = db_4vector[-1] / db_4vector[0]
        lorentz_boost = LorentzBoost(beta_z)

        # Align the momentum of W in the +z direction using two rotations
        boosted_wboson_4v = lorentz_boost.apply_transformation(wboson_4v)

        # 1. rotation around z axis by an angle -theta
        theta_angle = -np.arctan2(boosted_wboson_4v[2], boosted_wboson_4v[1])
        rotation_zaxis = LorentzRotationZ(theta_angle)

        # Leave the vector only in the z and x plane
        rotated_wboson_4v = rotation_zaxis.apply_transformation(boosted_wboson_4v)

        # 2. rotation around the y axis by an angle phi
        alpha = -np.arctan2(rotated_wboson_4v[1], rotated_wboson_4v[3])  # -atan2(px, pz)
        rotation_yaxis = LorentzRotationY(alpha)

        # Leave the boost direction in the +x and z plane (rotation around z by pi if needed)
        boost_direction = np.array([0, 0, 0, 1 if beta_z > 0 else -1])
        boost_direction = rotation_yaxis.apply_transformation(boost_direction)
        theta_last = -np.arctan2(boost_direction[2], boost_direction[1])  # -atan2(by, bx)
        rotation_boost = LorentzRotationZ(theta_last)

        # All Lorentz transformations that must be applied (order must be respected)
        final_transf = CompositeLorentzTransformations([lorentz_boost, rotation_zaxis, rotation_yaxis, rotation_boost])

        # Applies to all particles
        particles.apply_lorentz_transf(final_transf)

        return particles

