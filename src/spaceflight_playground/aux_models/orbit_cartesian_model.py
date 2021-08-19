#!/usr/bin/python3

import numpy as np
import casadi as cas

##
# @class orbit_cartesian_model
# @brief Model of a 2-dimensional point mass spacecraft in cartesian coordinates.
##
class orbit_cartesian_model:

    ##
    # @brief Initialization procedure
    ##
    def __init__(self):
        # Universe parameters

        # Gravitational constant (m^3/(kg*s^2))
        self.G = 6.67408 * 10**-11
        # Moon mass (kg)
        self.M = 7.348 * 10**22
        # Std. gravitational parameter (m^3/s^2)
        self.mu = self.G * self.M
        # Moon radius (m)
        self.R = 1737.5 * 10**3

        # Spacecraft parameters
        
        # Initial mass (kg)
        self.m0 = 20.0 * 10**3
        # Empty mass (kg)
        self.me = 1.0 * 10**3
        # Maximum thrust (N)
        self.u_max = 30.0 * 10**4
        # Fuel consumption coefficient (kg/s)
        self.km = 100.0

        # Dynamics parameters

        # State and control dimension
        self.nx = 5
        self.nu = 2
        # Scaling factor
        self.scale = cas.vertcat(
            1e-3, # m -> km
            1e-3, # m -> km
            1e-3, # m/s -> km/s
            1e-3, # m/s -> km/s
            1     # kg -> kg
        )
        # Un-scaling factor
        self.unscale = self.scale**-1
        # Initial state
        self.x0 = cas.vertcat(
            self.R,
            0.0,
            0.0,
            0.0,
            self.m0
        )
        # Scaled initial state
        self.x0_scaled = self.x0 * self.scale


    ##
    # @brief The ODE xdot = f(x,u) of the spacecraft
    # @param x The state (cartesian coordinates):
    #          - x-coordinate (m)
    #          - y-coordinate (m)
    #          - x-velocity (m/s)
    #          - y-velocity (m/s)
    #          - mass (kg)
    # @param u The controls
    #          - Thrust in x-direction (N)
    #          - Thrust in y-direction (N)
    # @return The derivative of the state
    ##
    def dynamics(self, x, u):

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
        xAcc = (T_x * self.u_max / m) - (xPos * coeff_grav)
        yAcc = (T_y * self.u_max / m) - (yPos * coeff_grav)

        # Compute mass derivative
        mDot = - self.km * (T_x**2 + T_y**2)
        
        # Stack the derivatives
        xdot = cas.vertcat(
            xVel,
            yVel,
            xAcc,
            yAcc,
            mDot
        )

        return xdot

    ##
    # @brief The scaled ODE of the spacecraft
    # @param x The scaled state (cartesian coordinates):
    #          - x-coordinate (km)
    #          - y-coordinate (km)
    #          - x-velocity (km/s)
    #          - y-velocity (km/s)
    #          - mass (kg)
    # @param u The controls
    #          - Thrust in x-direction (N)
    #          - Thrust in y-direction (N)
    # @return The derivative of the state
    ##
    def dynamics_scaled(self, x_scaled, u):
        # Unscale back to (m, m/s)
        x = x_scaled * self.unscale

        # Feed the ODE with the downscaled state
        xdot = self.dynamics(x, u)

        # Scale back to (km, km/s) and return
        xdot_scaled = xdot * self.scale
        return xdot_scaled


##
# Execute this script to run the model and to generate a trajectory
##
if __name__ == '__main__':

    import sys, os
    sys.path.append(os.path.realpath('../../../'))
    sys.path.append(os.getcwd())

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
    print("Initial mass: " + str(spacecraft.m0) + " kg")
    print("Empty mass: " + str(spacecraft.me) + " kg")
    print("Max. thrust: " + str(spacecraft.u_max) + " N")
    print("Fuel consumption coeff:" + str(spacecraft.km) + " kg/s")
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