import numpy as np
from typing import Tuple
from spaceflight_playground.constants import Universe
from spaceflight_playground.conversion import pol2cart, state_pol2cart


class KeplerOrbit:
    def __init__(self):
        """Representation of a keplerian orbit. Convenience class
        for easy computation of kepler elements from cartesian state vector.

        Coordinate system definitions:
        - X axis is horizontal, points right (vernal equinox)
        - Y axis is vertical, points upwards
        - Z axis points out of the screen (Body's rotational axis)
        """
        # Universe parameters
        # Std. gravitational parameter (m^3/s^2)
        self.mu = Universe.gravitational_constant * Universe.moon_mass

        # Orbital elements
        # Eccentricity vector (in 2d only the first two elements are non-zero)
        self.eccentricity = np.array([0, 0, 0])
        # Specific orbital angular momentum (in 2d only the last element is non-zero)
        self.specific_angular_momentum = np.array([0, 0, 0])

    def __str__(self):
        return f"eccentricity = {self.eccentricity}, angular momentum = {self.specific_angular_momentum}"

    def setOrbitalElements(self, eccentricity: np.ndarray, specific_angular_momentum: float) -> None:
        """Setter for the orbital elements
        :param eccentricity: The orbit's eccentricity.
        :param specific_angular_momentum: The orbit's specific orbital angular momentum.
        :return: None
        """
        self.eccentricity = eccentricity
        self.specific_angular_momentum = specific_angular_momentum


    def fromEllipseParams(
            self,
            eccentricity: float,
            angle: float,
            big_semi_major_axis: float
    ) -> Tuple[float, float]:
        """Computes kepler elements from common ellipse parameters.
        :param eccentricity: Eccentricity (Scalar value)
        :param angle: Angle of the orbit
        :param big_semi_major_axis: Big semi-major axis
        :return: Eccentricity vector, specific orbital angular momentum
        """
        # Compute the norm of specific orbital angular momentum
        semi_latus_rectum = big_semi_major_axis * (1 - eccentricity ** 2)
        specific_angular_momentum = np.sqrt(semi_latus_rectum * self.mu)

        # Compute and update kepler elements in 2d, lifted into 3d
        self.specific_angular_momentum = np.array([0, 0, specific_angular_momentum])
        self.eccentricity = np.array([eccentricity * np.cos(angle), eccentricity * np.sin(angle), 0])
        return self.eccentricity, self.specific_angular_momentum


    def fromCartesianState (
            self,
            position_x: float,
            position_y: float,
            velocity_x: float,
            velocity_y: float
    ) -> Tuple[float, float]:
        """Computes kepler elements from a cartesian state vector.
        Updates the internal orbit parameters.
        :param position_x: The x-position of the orbitting object.
        :param position_y: The y-position of the orbitting object.
        :param velocity_x: The x-velocity of the orbitting object.
        :param velocity_y: The y-velocity of the orbitting object.
        :return: Eccentricity vector, specific orbital angular momentum.
        """
        # Expand dimensionality and norm the position
        position = np.array([position_x, position_y, 0])
        velocity = np.array([velocity_x, velocity_y, 0])
        normed_position = position / np.linalg.norm(position)

        # Compute and update kepler elements
        self.specific_angular_momentum = np.cross(position, velocity)
        self.eccentricity = np.cross(velocity, self.specific_angular_momentum) / self.mu - normed_position
        return self.eccentricity, self.specific_angular_momentum


    def fromPolarState (
            self,
            distance: float,
            angle: float,
            radial_velocity: float,
            angular_velocity: float
    ) -> Tuple[float, float]:
        """Computes kepler elements from polar state vector.
        Updates the internal orbit parameters.
        :param distance: The orbiter's distance to the origin.
        :param angle: The orbiter's angle, measured counter-clockwise from the positive x-axis.
        :param radial_velocity: The orbiter's velocity away from the origin.
        :param angular_velocity: The orbiter's angular velocity.
        :return: Eccentricity vector, specific orbital angular momentum.
        """
        # Convert polar coordinates to cartesian coordinates and compute kepler elements
        state_polar = np.array([distance, angle, radial_velocity, angular_velocity])
        state_cartesian = state_pol2cart(state_polar)
        return self.fromCartesianState(state_cartesian[0], state_cartesian[1], state_cartesian[2], state_cartesian[3])


    def discretize (self, num_samples=360) -> np.ndarray:
        """Discretizes the orbit into a bunch of positions.
        :param num_samples: Number of discrete samples.
        :return: A Nx2 vector representing the discretized orbit.
        """
        eccentricity_norm = np.linalg.norm(self.eccentricity)
        angular_momentum_norm = np.linalg.norm(self.specific_angular_momentum)
        semi_latus_rectum = angular_momentum_norm**2 / self.mu

        # Compute rotation of the ellipse
        ellise_angle = np.arctan2(self.eccentricity[1], self.eccentricity[0])
        
        # Prepare samples container
        samples = np.zeros((num_samples, 2))

        # Discretize the ellipse
        angles = np.linspace(0, 2 * np.pi, num_samples)
        for k in range(num_samples):
            # Get current angle
            angle = angles[k]
            
            # Compute radius at current angle ("orbit equation")
            distance = semi_latus_rectum / (1 + eccentricity_norm * np.cos(angle))

            # Rotate the orbit clockwise
            angle = angle + ellise_angle

            # Convert into cartesian coordinates
            position = pol2cart(angle, distance)

            # Store coordinates in container
            samples[k, 0] = position[0]
            samples[k, 1] = position[1]

        # Return the discretized orbit
        return samples
