#!/usr/bin/python3

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


##
# Execute this script to run the model and to generate a trajectory
##    
if __name__ == '__main__':

    # Import plotting library and runge kutta 4 integrator    
    import matplotlib.pyplot as plt
    from src.spaceflight_playground.rk4step import rk4step_ode
    from src.spaceflight_playground.aux_models.kepler_orbit import KeplerOrbit as orbit
    from src.spaceflight_playground.conversion import *
    from src.spaceflight_playground.xml_writer import *

    # Create a spacecraft instance
    spacecraft = orbit_polar_model()

    # Print some parameters
    print("Universe parameters:")
    print("Grav. const. G: " + str(spacecraft.gravitational_constant) + " m^3/(kg*s^2)")
    print("Moon mass M: " + str(spacecraft.moon_mass) + " kg")
    print("Moon radius R: " + str(spacecraft.moon_radius) + " m")
    
    print("Spacecraft parameters:")
    print("Initial mass: " + str(spacecraft.full_mass) + " kg")
    print("Empty mass: " + str(spacecraft.empty_mass) + " kg")
    print("Max. thrust: " + str(spacecraft.max_thrust) + " N")
    print("Fuel consumption coeff:" + str(spacecraft.km) + " kg/s")
    print("Initial state: " + str(spacecraft.x0))
    print("Initial state (scaled): " + str(spacecraft.x0_scaled))
    print("State scale vector: " + str(spacecraft.scale))
    print("State unscale vector: " + str(spacecraft.unscale))
    print("Number of states: " + str(spacecraft.num_states))
    print("Number of controls: " + str(spacecraft.num_forces))

    # Simulation parameters
    T = 600.0
    N = 100
    DT = T/N

    # Integration parameters 
    nn = 10     # Integrator steps per step
    h = DT/nn   # Step size of integrator step

    # Create system model with CasADi
    x = cas.MX.sym('x', spacecraft.num_states, 1)
    u = cas.MX.sym('u', spacecraft.num_forces, 1)
    f = cas.Function('f', [x,u], [spacecraft.ode_scaled(x, u)])

    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        Xk = rk4step_ode(f, Xk, u, h)

    F = cas.Function('F', [x,u], [Xk], ['x','u'], ['xk'])

    # Choose controls for simulation
    us = np.zeros((spacecraft.num_forces, N))
    n_r_stop = 60
    n_theta_stop = 85
    
    us_r = np.zeros(N)
    us_r[0:n_r_stop] = 0.1075 * np.ones(n_r_stop)
    us[0,:] = us_r
    
    us_theta = np.zeros(N)
    us_theta[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
    us[1,:] = us_theta
    
    # Simulate the system
    xs = cas.DM.zeros((spacecraft.num_states, N + 1))
    xs[:,0] = spacecraft.x0_scaled

    for k in range(1,N+1):
        xs[:,k] = F(xs[:,k-1],us[:,k-1])


    # Unscale trajectory
    for k in range(N+1):
        xs[:,k] = xs[:,k] * spacecraft.unscale


    xs = xs.full()

    # Prepare plotting
    tAxis = np.linspace(0, T-DT, N+1)
    plt.figure(1)

    # Reference altitude [km] & angular velocity [µrad]
    altitude_ref = 20
    angVel_ref = 1e3 * np.sqrt(spacecraft.mu / (spacecraft.moon_radius + 1e3 * altitude_ref) ** 3)

    print("Reference altitude: " + str(altitude_ref))
    print("Reference angular velocity: " + str(angVel_ref))

    # Plot 
    plt.subplot(321)
    plt.plot(tAxis, xs[0,:])
    plt.plot(tAxis, altitude_ref * np.ones((tAxis.size,1)))
    plt.ylabel('Altitude [km]')

    plt.subplot(322)
    plt.plot(tAxis, xs[1,:])
    plt.plot(tAxis, spacecraft.x0[1].full() * np.ones((tAxis.size,1)))
    plt.ylabel('Angle [murad]')

    plt.subplot(323)
    plt.plot(tAxis, xs[2,:])
    plt.plot(tAxis, spacecraft.x0[2].full() * np.ones((tAxis.size,1)))
    plt.ylabel('Radial vel. [km/s]')

    plt.subplot(324)
    plt.plot(tAxis, xs[3,:])
    plt.plot(tAxis, angVel_ref * np.ones((tAxis.size,1)))
    plt.ylabel('Angular vel. [murad/s]')

    plt.subplot(325)
    plt.plot(tAxis, xs[4,:])
    plt.ylabel('Mass [kg]')

    plt.show()

    # Convert trajectory to cartesian coordinates
    for k in range(N+1):
        xs[0,k] = xs[0,k] + spacecraft.moon_radius
    
    # Unscale trajectory
    #for k in range(N+1):
    #    xs[:,k] = xs[:,k] * spacecraft.unscale


    xs_cart, us_cart = traj_pol2cart(xs[0:4,:], us)

    print(xs_cart[0,:])
    print(xs_cart[1,:])

    # Create an orbit instance
    orb = orbit()
    orb.fromPolarState(
        distance=spacecraft.moon_radius + spacecraft.unscale[0] * altitude_ref,
        angle= 0.0,
        radial_velocity= 0.0,
        angular_velocity=spacecraft.unscale[3] * angVel_ref
    )

    # Write trajectory to xml file
    write_to_xml(
        filename = 'orbit_animator_trajectory.xml',
        T = T,
        N = N,
        e =orb.eccentricity[0:2],
        h = orb.specific_angular_momentum[2],
        xPoses = xs_cart[0,:],
        yPoses = xs_cart[1,:],
        xVelos = xs_cart[2,:],
        yVelos = xs_cart[3,:],
        masses = xs[4,:],
        xForces = us_cart[0,:],
        yForces = us_cart[1,:]
    )