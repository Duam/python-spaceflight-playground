#!/usr/bin/python3

import casadi as cas
import numpy as np

class 2d_liftoff_model:
    ## 
    # @brief Initialization procedure
    ##
    def __init__(self):
        # Distance from base to COM (m)
        self.L = 20
        # Distance from base to COP (m)
        self.l = 10
        # Maximum thrust (N)
        self.maxThrust = 100
        # Maximum gimbal angle (rad)
        self.maxGimbal = np.pi / 2
        # Gravitational constant (m/s^2)
        self.g = 9.81
        # Total mass (kg)
        self.m = 1000
        # Inertia tensor (kg*m^2)
        self.I = 1
    
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
    # @param d (disturbance)
    #    Wind force (N)
    # @return State derivative
    ##
    def dynamics(self, x, u, d):
        # Get states
        xPos = x[0]
        yPos = x[1]
        xVel = x[2]
        yVel = x[3]
        ang = x[4]
        angVel = x[5]

        # Get controls
        T = u[0]
        mu = u[1]

        # Get disturbance (only in horizontal direction)
        F_W = d

        # Compute some intermediate values
        sa = cas.sin(ang)
        ca = cas.cos(ang)
        sm = cas.sin(mu)
        cm = cas.cos(mu)

        # Compute the thruster forces
        F_T = T * cas.vertcat(
            sm * sa + cm * ca, 
            sm * ca - cm * sa
        )

        # Compute the gravitational force
        F_G = cas.vertcat(
            0,
            - self.m * self.g
        )

        # Compute the torque around base
        Torque = self.L * F_W * ca + self.l * F_G * sa

        # Compute the linear accelerations
        xAcc = (F_T[0] + F_W) / self.m
        yAcc = (F_T[1] + F_G) / self.m

        # Compute the angular acceleration
        angAcc = Torque / self.I

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