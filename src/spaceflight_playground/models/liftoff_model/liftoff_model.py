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
        self.g = self.params['g']

        # Spacecraft parameters

        # Distance from COM to base (m)
        self.L = self.params['L']
        # Distance from COM to COP (m)
        self.l = self.params['l']
        # Maximum thrust (N)
        self.maxThrust = self.params['maxThrust']
        # Maximum gimbal angle in both directions (rad)
        self.maxGimbal = self.params['maxGimbal']
        # Total mass (kg)
        self.m = self.params['m']
        # Inertia tensor (kg*m^2)
        self.I = self.params['I']

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

