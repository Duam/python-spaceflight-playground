#!/usr/bin/python3

import numpy as np
import casadi as cas

##
# @class orbit_2d_polar_model
# @brief TODO
##
class orbit_2d_polar_model:

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
        self.km = 10**2
        # Initial state
        self.x0 = cas.vertcat(
            0.0,
            np.pi / 2.0,
            0.0,
            0.0,
            self.m0
        )

        # Dynamics parameters
        self.nx = 5
        self.nu = 2

        # Scaling factor
        # km -> m
        # nrad -> rad
        self.scaledown = cas.vertcat(
            1e3,
            1e-6,
            1e3,
            1e-6,
            1
        )
        self.scaleup = self.scaledown**-1 # Does this work?

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
        acc_theta_thrust = theta * self.u_max / (m * (r+self.R))
        # Coriolis term 
        acc_theta_corio = - 2 * rDot * thetaDot / (r+self.R)
        # Unify terms
        thetaDDot = acc_theta_thrust + acc_theta_corio

        # Compute mass derivative
        mDot = - self.km * (T_r**2 + T_theta**2)
        
        # Stack the derivatives
        return cas.vertcat(
            rDot,
            thetaDot,
            rDDot,
            thetaDDot,
            mDot
        )

    ##
    # @brief The scaled ODE of the spacecraft. 
    # @param x The scaled state (polar coordinates):
    #          - Altitude (km)
    #          - Angle from horizontal axis (nrad)
    #          - Radial velocity (km/s)
    #          - Angular velocity (nrad/s)
    #          - mass (kg)
    # @param u The controls
    #          - Thrust in radial direction (N)
    #          - Thrust in tangential direction (N)
    # @return The scaled derivative of the state
    ##
    ##
    def dynamics_scaled(self, x_upscaled, u):
        # Downscale state from (km,nrad,..) to (m,rad,..)
        x_downscaled = x_upscaled * self.scaledown

        # Feed the ODE with the downscaled state
        xdot_downscaled = self.dynamics(x_downscaled, u)

        # Upscale the state derivative
        return xdot_downscaled * self.scaleup


##
# Execute this script to test the model
##    
if __name__ == '__main__':

    import sys, os
    sys.path.append(os.path.realpath('../../'))
    sys.path.append(os.getcwd())

    # Import plotting library and runge kutta 4 integrator    
    import matplotlib.pyplot as plt
    from integrators.rk4step import rk4step_ode

    # Create a spacecraft instance
    spacecraft = orbit_2d_polar_model()

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
    print("State upscale vector: " + str(spacecraft.scaleup))
    print("State downscale vector: " + str(spacecraft.scaledown))
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
    ode_scl = cas.Function('ode_scl', [x,u], [spacecraft.dynamics_scaled(x,u)])

    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        #Xk = rk4step_ode(ode, Xk, u, h)
        Xk = rk4step_ode(ode_scl, Xk, u, h)

    #ode_d = cas.Function('ode_d', [x,u], [Xk], ['x','u'], ['xk'])
    ode_scl_d = cas.Function('ode_scl_d', [x,u], [Xk], ['x','u'], ['xk'])

    # Choose controls for simulation
    us = np.zeros((N,spacecraft.nu))
    n_r_stop = 60
    n_theta_stop = 85
    
    us_r = np.zeros(N)
    us_r[0:n_r_stop] = 0.1705 * np.ones(n_r_stop)
    us[:,0] = us_r
    
    us_theta = np.zeros(N)
    us_theta[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
    us[:,1] = us_theta
    
    # Simulate the system
    xs = cas.DM.zeros((N, spacecraft.nx))
    xs[0,:] = spacecraft.x0

    for k in range(1,N):
        xs[k,:] = ode_scl_d(xs[k-1,:],us[k-1,:])

    xs = xs.full()

    # Prepare plotting
    tAxis = np.linspace(0, T-DT, N)
    plt.figure(1)

    # Reference altitude [km] & angular velocity [µrad]
    altitude_ref = 20
    angVel_ref = 10**6 * np.sqrt(spacecraft.mu / (spacecraft.R + 10**3 * altitude_ref)**3)

    # Plot 
    plt.subplot(321)
    plt.plot(tAxis, xs[:,0])
    plt.plot(tAxis, altitude_ref * np.ones((tAxis.size,1)))
    plt.ylabel('Altitude [km]')

    plt.subplot(322)
    plt.plot(tAxis, xs[:,1])
    plt.plot(tAxis, spacecraft.x0[1].full() * np.ones((tAxis.size,1)))
    plt.ylabel('Angle [murad]')

    plt.subplot(323)
    plt.plot(tAxis, xs[:,2])
    plt.plot(tAxis, spacecraft.x0[2].full() * np.ones((tAxis.size,1)))
    plt.ylabel('Radial vel. [km/s]')

    plt.subplot(324)
    plt.plot(tAxis, xs[:,3])
    plt.plot(tAxis, angVel_ref * np.ones((tAxis.size,1)))
    plt.ylabel('Angular vel. [murad/s]')

    plt.subplot(325)
    plt.plot(tAxis, xs[:,4])
    plt.ylabel('Mass')

    plt.show()