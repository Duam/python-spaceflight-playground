import numpy as np
from dataclasses import dataclass
from typing import Tuple, List

@dataclass
class PolarState:
    distance: float  # The distance to the origin [m].
    angle: float  # The angle to the positive x-axis, measured counter-clockwise [rad].
    radial_velocity: float  # The velocity component in radial direction [m/s].
    angular_velocity: float  # The velocity component in tangential direction [rad/s].
    def as_numpy_vector(self):
        return np.array([distance, angle, radial_velocity, angular_velocity])

@dataclass
class CartesianState:
    x_position: float
    y_position: float
    x_velocity: float
    y_velocity: float
    def as_numpy_vector(self):
        return np.array([x_position, y_position, x_velocity, y_velocity])

@dataclass
class PolarForce:
    radial: float
    angular: float
    def as_numpy_vector(self):
        return np.array([radial, angular])

@dataclass
class CartesianForce:
    x: float
    y: float
    def as_numpy_vector(self):
        return np.array([x, y])


def cart2pol (x_position: float, y_position: float) -> Tuple[float, float]:
    """Converts cartesian coordinates to polar coordinates.
    :param x_position: The x coordinate in the cartesian frame.
    :param y_position: The y coordinate in the cartesian frame.
    :return: The polar coordinates [angle, distance].
    """
    angle = np.arctan2(y, x)
    distance = np.hypot(x, y)
    return angle, distance


def pol2cart (angle: float, distance: float) -> Tuple[float, float]:
    """Converts polar coordinates to cartesian coordinates.
    :param angle: The angle coordinate in the polar frame.
    :param distance: The distance coordinate in the polar frame.
    :return: The cartesian coordinates [x-position, y-position].
    """
    x_position = distance * np.cos(angle)
    y_position = distance * np.sin(angle)
    return x_position, y_position


def state_pol2cart(state_pol: PolarState) -> CartesianState:
    """Converts a state vector in polar coordinates to cartesian coordinates.
    :param state_pol: The state in polar coordinates.
    :return: The state in cartesian coordinates.
    """
    # Compute position in cartesian frame and expand into 3d (for cross-product)
    x_position, y_position = pol2cart(state_pol.angle, state_pol.distance)

    # Expand quantities into 3d (for cross-product)
    position = np.array([x_position, y_position, 0])
    angular_velocity = np.array([0, 0, state_pol.angular_velocity])

    # Compute radial and rotational velocities in cartesian frame
    radial_velocity = np.array([state_pol.radial_velocity, 0, 0])
    tangential_velocity = np.cross(angular_velocity, position)
    velocity = radial_velocity + tangential_velocity

    return CartesianState(position[0], position[1], velocity[0], velocity[1])


def force_pol2cart(angle: float, force_pol: PolarForce) -> CartesianForce:
    """Converts a force from polar coordinates to cartesian coordinates.
    :param angle: The angle of the corotating frame.
    :param force_pol: The force in polar coordinates.
    :return: The force in cartesian coordinates
    """
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    return CartesianForce(np.dot(rotation, force_pol.as_numpy_vector()))


def traj_pol2cart (
        states_pol: List[PolarState],
        forces_pol: List[PolarForce],
) -> Tuple[List[CartesianState], List[CartesianForce]]:
    """Converts a state and control trajectory from polar coordinates to cartesian coordinates.
    :param states_pol: The state trajectory in polar coordinates.
    :param controls_pol: The force trajectory in polar coordinates.
    :return: The state trajectory in cartesian coordinates, the force trajectory in cartesian coordinates.
    """
    # Check if trajectories are compatible
    assert states_pol.size == forces_pol.size + 1,\
        f"There has to be one more state than there" \
        f"are controls! Nx = {states_pol.size}, Nu = {forces_pol.size}"

    # Convert states and forces to cartesian coordinates
    states_cart = [state_pol2cart(state_pol) for state_pol in states_pol]
    forces_cart = [force_pol2cart(states_pol[k].angle, forces_pol[k]) for k in range(forces_pol.size)]
    return states_cart, forces_cart
