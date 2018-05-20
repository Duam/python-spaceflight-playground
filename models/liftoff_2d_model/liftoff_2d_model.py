#!/usr/bin/python3

import casadi as cas
import numpy as np

## 
# @class liftoff_2d_model
# @brief Model of a 2d rocket during liftoff while 
#        considering very simple aerodynamics.
##
class liftoff_2d_model:
  
    ## 
    # @brief Initialization procedure
    # @param L Distance from base to COM (m)
    # @param l Distance from base to COP (m)
    # @param maxThrust Maximum thrust (N)
    # @param maxGimbal Maximum thrust gimbal angle (rad)
    # @param m Total mass (kg)
    # @param I Inertia tensor (around z-axis) (kg*m^2)
    ##
    def __init__(self, L=20, l=10, maxThrust=100, maxGimbal=np.pi/2, m=1000, I=1):
        # Gravitational constant (m/s^2)
        self.g = 9.81

        # Distance from base to COM (m)
        self.L = L
        # Distance from base to COP (m)
        self.l = l
        # Maximum thrust (N)
        self.maxThrust = maxThrust
        # Maximum gimbal angle (rad)
        self.maxGimbal = maxGimbal
        # Total mass (kg)
        self.m = m
        # Inertia tensor (kg*m^2)
        self.I = I
    
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

        # Get disturbing wind force (only in horizontal direction)
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