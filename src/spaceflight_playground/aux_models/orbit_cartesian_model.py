import numpy as np
import casadi as cas
from dataclasses import dataclass
from spaceflight_playground.constants import Universe
from typing import Union


class CartesianOrbiter:
    def __init__(self):
        """Model of a 2-dimensional point mass spacecraft in cartesian coordinates.

        Assumptions/Simplifications:
        - No atmosphere, because it's the moon
        - The moon doesn't rotate
        - Spacecraft has no rotation
        - Thrusters can fire in any direction
        - The moon is a point mass and perfectly circular
        """
        # Standard gravitational parameter [m^3/s^2]
        self.mu = Universe.gravitational_constant * Universe.moon_mass
        self.full_mass = 20.0 * 1e3  # [kg]
        self.empty_mass = 1.0 * 1e3  # [kg]
        self.max_thrust = 30.0 * 1e4  # [N]
        self.fuel_consumption_per_second = 100.0  # [kg/s]

        # Scaling factor
        self.scale = cas.vertcat(
            1e-3,  # m -> km
            1e-3,  # m -> km
            1e-3,  # m/s -> km/s
            1e-3,  # m/s -> km/s
            1      # kg -> kg
        )
        # Un-scaling factor
        self.unscale = self.scale**-1

        # Initial state
        self.x0 = cas.vertcat(
            Universe.moon_radius,
            0.0,
            0.0,
            0.0,
            self.full_mass
        )
        # Scaled initial state
        self.x0_scaled = self.x0 * self.scale


    def dynamics(self, state, thrust):
        """The ODE describing the dynamics of the spacecraft.
        State:
        - x-coordinate [m]
        - y-coordinate [m]
        - x-velocity [m/s]
        - y-velocity [m/s]
        - mass [kg]
        Thrust:
        - Force in x-direction [N]
        - Force in y-direction [N]
        :param state: The current state.
        :param thrust: The currently acting force.
        :return: The current state derivative.
        """

        # Get states
        xPos = x[0] # x-Position (m)
        yPos = x[1] # y-Position (m)
        xVel = x[2] # x-velocity (m/s)
        yVel = x[3] # y-velocity (m/s)
        m = x[4]    # Mass (kg)

        # Get controls
        T_x = u[0] # x-thrust (N)
        T_y = u[1] # y-thrust (N)

        # Precompute gravitational force coefficient
        coeff_grav = self.mu / cas.sqrt(xPos**2 + yPos**2)**3

        # Compute accelerations
        xAcc = (T_x * self.max_thrust / m) - (xPos * coeff_grav)
        yAcc = (T_y * self.max_thrust / m) - (yPos * coeff_grav)

        # Compute mass derivative
        mDot = - self.fuel_consumption_per_second * (T_x ** 2 + T_y ** 2)
        
        # Stack the derivatives
        xdot = cas.vertcat(
            xVel,
            yVel,
            xAcc,
            yAcc,
            mDot
        )

        return xdot


    def dynamics_scaled(self, state_scaled, thrust):
        """The scaled ODE of the spacecraft.
        :param state_scaled: The scaled state.
        :param thrust: The thrust vector applied to the spacecraft.
        :return: The scaled state derivative.
        """
        state = state_scaled * self.unscale
        state_derivative = self.dynamics(state, thrust)
        return state_derivative * self.scale


##
# Execute this script to run the model and to generate a trajectory
##
if __name__ == '__main__':

    # Import plotting library and runge kutta 4 integrator    
    import matplotlib.pyplot as plt
    from src.spaceflight_playground.rk4step import rk4step_ode

    # Create a spacecraft instance
    spacecraft = orbit_cartesian_model()

    # Print some parameters
    print("Universe parameters:")
    print("Grav. const. G: " + str(spacecraft.G) + " m^3/(kg*s^2)")
    print("Moon mass M: " + str(spacecraft.M) + " kg")
    print("Moon radius R: " + str(spacecraft.R) + " m")
    
    print("Spacecraft parameters:")
    print("Initial mass: " + str(spacecraft.full_mass) + " kg")
    print("Empty mass: " + str(spacecraft.empty_mass) + " kg")
    print("Max. thrust: " + str(spacecraft.max_thrust) + " N")
    print("Fuel consumption coeff:" + str(spacecraft.fuel_consumption_per_second) + " kg/s")
    print("Initial state: " + str(spacecraft.x0))
    print("Initial state (scaled): " + str(spacecraft.x0_scaled))
    print("State scale vector: " + str(spacecraft.scale))
    print("State unscale vector: " + str(spacecraft.unscale))
    print("Number of states: " + str(spacecraft.nx))
    print("Number of controls: " + str(spacecraft.nu))

    # Simulation parameters
    T = 600.0
    N = 100
    DT = T/N

    # Integration parameters 
    nn = 10     # Integrator steps per step
    h = DT/nn   # Step size of integrator step

    # Create system model with CasADi
    x = cas.MX.sym('x', spacecraft.nx, 1)
    u = cas.MX.sym('u', spacecraft.nu, 1)
    #ode = cas.Function('ode', [x,u], [spacecraft.dynamics(x,u)])
    f = cas.Function('f', [x,u], [spacecraft.dynamics_scaled(x,u)])

    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        Xk = rk4step_ode(f, Xk, u, h)
    F = cas.Function('F', [x,u], [Xk], ['x','u'], ['xk'])

    # Choose controls for simulation
    us = np.zeros((N,spacecraft.nu))
    n_x_stop = 25
    n_y_stop = 50
    
    us[0:n_x_stop,0] = 0.11 * np.ones(n_x_stop)
    us[0:n_y_stop,1] = 0.5 * np.ones(n_y_stop)

    # Simulate the system (x = state vector)
    xs = cas.DM.zeros((N, spacecraft.nx))
    xs[0,:] = spacecraft.x0_scaled

    for k in range(1,N):
        xs[k,:] = F(xs[k-1,:],us[k-1,:])

    xs = xs.full()

    # Prepare plotting
    tAxis = np.linspace(0, T-DT, N)
    plt.figure(1)

    # Plot 
    plt.subplot(311)
    plt.plot(xs[:,0], xs[:,1])
    plt.xlabel('x-position [km]')
    plt.ylabel('y-position [km]')

    plt.subplot(312)
    plt.plot(xs[:,2], xs[:,3])
    plt.xlabel('x-velocity [km/s]')
    plt.ylabel('y-velocity [km/s]')

    plt.subplot(313)
    plt.plot(tAxis, xs[:,4])
    plt.ylabel('Mass [kg]')

    plt.show()