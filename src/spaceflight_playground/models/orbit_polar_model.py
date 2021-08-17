#!/usr/bin/python3

import casadi as cas

"""
# Model
This model is that of a two dimensional point mass rocket. Its thrusters can fire in any direction.

## Inputs
- T_r: Scaled thrust in radial direction (Between 0 and 1)
- T_theta: Scaled thrust in tangential direction (Between 0 and 1)

## Disturbances
No disturbances

## States
- r: Altitude (km)
- theta: Angle from horizontal axis (microRad)
- rDot: Radial velocity (km per second)
- thetaDot: Angular velocity (microRad per second)
- m: Mass (kg)

# Assumptions/Simplifications:
- No atmosphere, because it's the moon
- The moon doesn't rotate
- Spacecraft has no rotation
- Thrusters can fire in any direction
- The moon is a point mass and perfectly circular
"""

##
# @class orbit_polar_model
# @brief Model of a 2-dimensional point mass spacecraft in polar coordinates.
##
class orbit_polar_model:

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
            1e+3, # rad -> millirad
            1e-3, # m/s -> km/s
            1e+3, # rad/s -> millirad/s
            1e-3  # kg -> t
        )
        self.unscale = self.scale**-1

        # Initial state
        self.x0 = cas.vertcat(
            0.0,
            0.0,
            0.0,
            0.0,
            self.m0
        )
        # Scaled initial state
        self.x0_scaled = self.x0 * self.scale

    ##
    # @brief The ODE xdot = f(x,u) of the spacecraft
    # @param x The state (polar coordinates):
    #          - Altitude (m)
    #          - Angle from horizontal axis (rad)
    #          - Radial velocity (m/s)
    #          - Angular velocity (rad/s)
    #          - mass (kg)
    # @param u The controls
    #          - Thrust in radial direction (N)
    #          - Thrust in tangential direction (N)
    # @return The derivative of the state
    ##
    def dynamics(self, x, u):
        
        # Get states
        r = x[0]        # Position (Altitude)
        theta = x[1]    # Position (angular)
        rDot = x[2]     # Velocity (radial)
        thetaDot = x[3] # Velocity (angular)
        m = x[4]        # Mass

        # Get controls
        T_r = u[0]     # Thrust (radial)
        T_theta = u[1] # Thrust (angular)

        # Compute acceleration in radial direction
        # Control term
        acc_r_thrust = T_r * self.u_max / m
        # Gravitaty term
        acc_r_grav = - self.mu / (r+self.R)**2
        # Centrifugal term
        acc_r_centri = (r+self.R) * thetaDot**2
        # Unify terms
        rDDot = acc_r_thrust + acc_r_grav + acc_r_centri

        # Compute angular acceleration (in tangential direction)
        # Control term
        acc_theta_thrust = T_theta * self.u_max / (m * (r+self.R))
        # Coriolis term 
        acc_theta_corio = - 2 * rDot * thetaDot / (r+self.R)
        # Unify terms
        thetaDDot = acc_theta_thrust + acc_theta_corio

        # Compute mass derivative
        mDot = - self.km * (T_r**2 + T_theta**2)

        # Stack the derivatives
        xdot = cas.vertcat(
            rDot,
            thetaDot,
            rDDot,
            thetaDDot,
            mDot
        )

        return xdot

    ##
    # @brief The scaled ODE of the spacecraft. 
    # @param x The scaled state (polar coordinates):
    #          - Altitude (km)
    #          - Angle from horizontal axis (microrad)
    #          - Radial velocity (km/s)
    #          - Angular velocity (microrad/s)
    #          - mass (kg)
    # @param u The controls
    #          - Thrust in radial direction (N)
    #          - Thrust in tangential direction (N)
    # @return The scaled derivative of the state
    ##
    ##
    def dynamics_scaled(self, x_scaled, u):
        # Unscale back to (m, rad, ..)
        x = x_scaled * self.unscale

        # Feed the ODE with the downscaled state
        xdot = self.dynamics(x, u)

        # Scale back to (km, microrad, ..) and return
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
    from src.spaceflight_playground.models.kepler_orbit import kepler_orbit as orbit
    from src.spaceflight_playground.conversion import *
    from src.spaceflight_playground.xml_writer import *

    # Create a spacecraft instance
    spacecraft = orbit_polar_model()

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
    f = cas.Function('f', [x,u], [spacecraft.dynamics_scaled(x,u)])

    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        Xk = rk4step_ode(f, Xk, u, h)

    F = cas.Function('F', [x,u], [Xk], ['x','u'], ['xk'])

    # Choose controls for simulation
    us = np.zeros((spacecraft.nu,N))
    n_r_stop = 60
    n_theta_stop = 85
    
    us_r = np.zeros(N)
    us_r[0:n_r_stop] = 0.1075 * np.ones(n_r_stop)
    us[0,:] = us_r
    
    us_theta = np.zeros(N)
    us_theta[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
    us[1,:] = us_theta
    
    # Simulate the system
    xs = cas.DM.zeros((spacecraft.nx,N+1))
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
    angVel_ref = 1e3 * np.sqrt(spacecraft.mu / (spacecraft.R + 1e3 * altitude_ref)**3)

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
        xs[0,k] = xs[0,k] + spacecraft.R
    
    # Unscale trajectory
    #for k in range(N+1):
    #    xs[:,k] = xs[:,k] * spacecraft.unscale


    xs_cart, us_cart = traj_pol2cart(xs[0:4,:], us)

    print(xs_cart[0,:])
    print(xs_cart[1,:])

    # Create an orbit instance
    orb = orbit()
    orb.fromPolarState(
        rho = spacecraft.R + spacecraft.unscale[0] * altitude_ref,
        theta = 0.0,
        rhoDot = 0.0,
        thetaDot = spacecraft.unscale[3] * angVel_ref
    )

    # Write trajectory to xml file
    write_to_xml(
        filename = 'trajectory.xml',
        T = T,
        N = N,
        e = orb.e[0:2],
        h = orb.h[2],
        xPoses = xs_cart[0,:],
        yPoses = xs_cart[1,:],
        xVelos = xs_cart[2,:],
        yVelos = xs_cart[3,:],
        masses = xs[4,:],
        xForces = us_cart[0,:],
        yForces = us_cart[1,:]
    )