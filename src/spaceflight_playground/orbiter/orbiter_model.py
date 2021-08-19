import casadi as cas
import numpy as np
from dataclasses import dataclass
from spaceflight_playground.conversion import PolarForce
from typing import Union


@dataclass
class PolarOrbiterState:
    vector: Union[cas.DM, cas.MX]
    def scale(self, factors: Union[cas.DM]):
        self.vector *= factors
        return self
    @property
    def distance(self):
        """The distance to the origin [m]."""
        return self.vector[0]
    @property
    def angle(self):
        """The angle to the positive x-axis, measured counter-clockwise [rad]."""
        return self.vector[1]
    @property
    def radial_velocity(self):
        """The velocity component in radial direction [m/s]."""
        return self.vector[2]
    @property
    def angular_velocity(self):
        """The velocity component in tangential direction [rad/s]."""
        return self.vector[3]
    @property
    def mass(self):
        """The total mass of the orbiter [kg]"""
        return self.vector[4]


@dataclass
class PolarOrbiterThrust:
    vector: Union[cas.MX]
    @property
    def radial(self):
        """The radial thrust component [N]."""
        return self.vector[0]
    @property
    def angular(self):
        """The angular thrust component [kg * rad/s^2]."""
        return self.vector[1]


@dataclass
class PolarOrbiterStateDerivative:
    vector: Union[cas.DM, cas.MX]
    def scale(self, factors: Union[cas.DM]):
        self.vector *= factors
        return self
    @property
    def radial_velocity(self):
        """The velocity component in radial direction [m/s]."""
        return self.vector[0]
    @property
    def angular_velocity(self):
        """The velocity component in tangential direction [rad/s]."""
        return self.vector[1]
    @property
    def radial_acceleration(self):
        """The acceleration component in radial direction [m/s^2]."""
        return self.vector[2]
    @property
    def angular_acceleration(self):
        """The acceleration component in tangential direction [rad/s^2]."""
        return self.vector[3]
    @property
    def mass(self):
        """The total mass of the orbiter [kg]"""
        return self.vector[4]


class PolarOrbiter:
    def __init__(self):
        """Model of a 2-dimensional point mass spacecraft in polar coordinates.

        Assumptions/Simplifications:
        - No atmosphere, because it's the moon
        - The moon doesn't rotate
        - Spacecraft has no rotation
        - Thrusters can fire in any direction
        - The moon is a point mass and perfectly circular
        """
        # Universe parameters
        self.gravitational_constant = 6.67408 * 1e-11  # [m^3/(kg*s^2)]
        self.moon_mass = 7.348 * 1e22  # [kg]
        self.moon_radius = 1737.5 * 1e3  # [m]

        # Spacecraft parameters
        self.full_mass = 20.0 * 1e3  # [kg]
        self.empty_mass = 1.0 * 1e3  # [kg]
        self.max_thrust = 30.0 * 1e4  # [N]
        self.fuel_consumption_per_second = 100.0  # [kg/s]

        # Scaling factor
        self.scale = cas.vertcat(
            1e-6, # m -> Mm
            1e+3, # rad -> millirad
            1e-3, # m/s -> km/s
            1e+6, # rad/s -> microrad/s
            1e-3  # kg -> t
        )
        self.unscale = self.scale**-1

    @property
    def num_states(self):
        return 5

    @property
    def num_forces(self):
        return 2



    def ode(self, state: PolarOrbiterState, force: PolarForce) -> PolarOrbiterStateDerivative:
        """The ODE describing the dynamics of the spacecraft.
        :param state: The current state in polar coordinates.
        :param force: The currently acting force in polar coordinates.
        :return: The state derivative.
        """
        # Get states and forces
        distance = state.distance
        radial_velocity = state.radial_velocity
        angular_velocity = state.angular_velocity
        mass = state.mass
        radial_thrust = force.radial
        angular_thrust = force.angular

        # Radial acceleration: Thrust term + gravity term + centrifugal term
        radial_acceleration = 0
        radial_acceleration += radial_thrust * self.max_thrust / mass
        radial_acceleration -= self.gravitational_constant * self.moon_mass / distance ** 2
        radial_acceleration += distance * angular_velocity ** 2

        # Angular acceleration: Thrust term + Coriolis term
        angular_acceleration = 0
        angular_acceleration += angular_thrust * self.max_thrust / (mass * distance)
        angular_acceleration -= 2 * radial_velocity * angular_velocity / distance

        # Mass flow
        mass_flow = - self.fuel_consumption_per_second * (radial_thrust**2 * angular_thrust**2)

        return PolarOrbiterStateDerivative(cas.vertcat(
            radial_velocity,
            angular_velocity,
            radial_acceleration,
            angular_acceleration,
            mass_flow,
        ))


    def ode_scaled(self, state_scaled: PolarOrbiterState, thrust: PolarForce) -> PolarOrbiterStateDerivative:
        """The scaled ODE of the spacecraft.
        :param state_scaled: The scaled state in polar coordinates.
        :param thrust: The thrust applied to the spacecraft in the corotating frame.
        :return: The scaled state derivative in polar coordinates.
        """
        state = state_scaled.scale(self.unscale)
        state_derivative = self.ode(state, thrust)
        return state_derivative.scale(self.scale)
