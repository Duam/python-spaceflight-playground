#!/usr/bin/python3

import casadi as cas
import numpy as np

## 
# @class liftoff_model
# @brief Model of a 2d rocket during liftoff while 
#        considering very simple aerodynamics.
##
class liftoff_model:
  
    ## 
    # @brief Initialization procedure
    ##
    def __init__(self, params = None):
        
        if (params == None):
            self.params = {
                'g': 9.81,
                'L': 20.0,
                'l': 10.0,
                'maxThrust': 300e3,
                'maxGimbal': np.pi / 8,
                'm': 20e3,
                'I': 1e4
            }

        # Universe parameters

        # Gravitational constant (m/s^2)
        self.g = params['g']

        # Spacecraft parameters

        # Distance from COM to base (m)
        self.L = params['L']
        # Distance from COM to COP (m)
        self.l = params['l']
        # Maximum thrust (N)
        self.maxThrust = params['maxThrust']
        # Maximum gimbal angle in both directions (rad)
        self.maxGimbal = params['maxGimbal']
        # Total mass (kg)
        self.m = params['m']
        # Inertia tensor (kg*m^2)
        self.I = params['I']

        # Dynamics parameters

        # Number of states, disturbances, controls
        self.nx = 6 
        self.nd = 1
        self.nu = 2
        # Initial state
        self.x0 = cas.vertcat(
            0.0, # x-position
            0.0, # y-position
            0.0, # x-velocity
            0.0, # y-velocity
            0.0, # angle (measured from vertical axis)
            0.0  # angular velocity
        )

        # Names of states, controls disturbances
        self.x_keys = [
            'xPos',
            'yPos',
            'xVel', 
            'yVel',
            'ang',
            'angVel'
        ]
        self.u_keys = [
            'thrust',
            'gimbal'
        ]
        self.d_keys = [
            'xForce'
        ]

    
    ##
    # @brief The model dynamics
    # @param x (state)
    #    x-position (m)
    #    y-position (m)
    #    x-velocity (m/s)
    #    y-velocity (m/s)
    #    angle (rad)
    #    angular velocity (rad/s)
    # @param u (input)
    #    Thrust (N)
    #    Gimbal angle (rad)
    #    Wind force disturbance (N)
    # @return State derivative
    ##
    def dynamics(self, x, u):
        # Get states
        xPos = x[0]
        yPos = x[1]
        xVel = x[2]
        yVel = x[3]
        ang = x[4]
        angVel = x[5]

        # Get controls
        thrust_percentage = u[0]
        thrust_angle = u[1]

        # Get disturbing wind force (only in horizontal direction)
        force_wind = u[2]

        # Compute some intermediate values
        thrust_magnitude = thrust_percentage * self.maxThrust
        sa = cas.sin(ang)
        ca = cas.cos(ang)
        st = cas.sin(thrust_angle)
        ct = cas.cos(thrust_angle)

        # == LINEAR FORCES ==
        force_gravity = - self.m * self.g
        force_thrust = thrust_magnitude * cas.vertcat(
            ct * sa - st * ca,
            ct * ca + st * sa
        )

        # Compute the linear accelerations
        xAcc = (force_thrust[0] + force_wind) / self.m
        yAcc = (force_thrust[1] + force_gravity) / self.m

        # == TORQUES ==
        torque_wind = self.l * ca * force_wind
        torque_thrust = self.L * st * thrust_magnitude
        
        # Compute the angular acceleration
        angAcc = (torque_thrust + torque_wind) / self.I

        # Stack the derivatives
        xdot = cas.vertcat(
            xVel,
            yVel,
            xAcc,
            yAcc,
            angVel,
            angAcc
        )

        return xdot


##
# Run this script to test the model and generate a trajectory
##
if __name__ == '__main__':

    import sys, os
    sys.path.append(os.path.realpath('../../'))
    sys.path.append(os.getcwd())

    # Import plotting library and runge kutta 4 integrator    
    import matplotlib.pyplot as plt
    from integrators.rk4step import rk4step_ode

    # Create a spacecraft instance
    spacecraft = liftoff_model()

    # Print some parameters
    print("Spacecraft parameters:")
    print("Grav. accel. g: " + str(spacecraft.g) + " m/s^2")
    print("Distance base to COM: " + str(spacecraft.L) + " m")
    print("Distance base to COP: " + str(spacecraft.l) + " m")
    print("Mass: " + str(spacecraft.m) + " kg")
    print("Max. thrust: " + str(spacecraft.maxThrust) + " N")
    print("Initial state: " + str(spacecraft.x0))
    print("Number of states: " + str(spacecraft.nx))
    print("Number of disturbances: " + str(spacecraft.nd))
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
    u = cas.MX.sym('u', spacecraft.nu + spacecraft.nd, 1)
    f = cas.Function('f', [x,u], [spacecraft.dynamics(x,u)])
    
    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        Xk = rk4step_ode(f, Xk, u, h)
    
    F = cas.Function('F', [x,u], [Xk], ['x','u'], ['xk'])

    # Choose controls for simulation
    wind_forces = np.zeros((N, spacecraft.nd))
    controls = np.zeros((N, spacecraft.nu))
    controls[:,0] = 0.7 * np.ones(N)
    controls[50:51,1] = 0.000001 * np.ones(1)

    us = np.append(controls, wind_forces, axis=1)

    # Test evaluation of the ODE
    print(us[0,:])
    xdot_test = spacecraft.dynamics(spacecraft.x0, us[0,:])
    print(xdot_test)

    # Simulate the system
    xs = cas.DM.zeros((N, spacecraft.nx))
    xs[0,:] = spacecraft.x0

    for k in range(1,N):
        xs[k,:] = F(xs[k-1,:], us[k-1,:])

    xs = xs.full()

    # Prepare plotting
    tAxis = np.linspace(0, T-DT, N)
    plt.figure(1)

    # Plot
    plt.subplot(321)
    plt.plot(tAxis, xs[:,0])
    plt.ylabel('x-pos [m]')

    plt.subplot(322)
    plt.plot(tAxis, xs[:,1])
    plt.ylabel('y-pos [m]')

    plt.subplot(323)
    plt.plot(tAxis, xs[:,2])
    plt.ylabel('x-vel [m/s]')

    plt.subplot(324)
    plt.plot(tAxis, xs[:,3])
    plt.ylabel('y-vel [m/s]')

    plt.subplot(325)
    plt.plot(tAxis, xs[:,4])
    plt.ylabel('angle [rad]')

    plt.subplot(326)
    plt.plot(tAxis, xs[:,5])
    plt.ylabel('angular vel. [rad/s]')

    plt.show()