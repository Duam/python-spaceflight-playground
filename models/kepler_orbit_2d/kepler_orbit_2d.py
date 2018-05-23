#!/usr/bin/python3

##
# @file kepler_orbit_2d.py
# @author Paul Daum
##

import numpy as np

##
# @class kepler_orbit_2d
# @brief Representation of a 2d keplerian orbit. Convenience class
#        for easy computation of kepler elements from cartesian state vector.
#        Should later be compatible with polar state vector.
##
class kepler_orbit_2d:

    ##
    # @brief Initialization procedure
    ##
    def __init__(self):
        # Universe parameters
        
        # Gravitational constant (m^3/(kg*s^2))
        self.G = 6.67408 * 10**-11
        # Moon mass (kg) TODO should be settable
        self.M = 7.348 * 10**22
        # Std. gravitational parameter (m^3/s^2)
        self.mu = self.G * self.M
        # Moon radius (m)
        self.R = 1737.5 * 10**3
        
        # Orbital elements

        # Eccentricity vector
        self.e = np.array([0,0])
        # Specific orbital angular momentum
        self.h = np.array([0,0])
        

    ##
    # @brief Computes kepler elements from cartesian state vector.
    #        Updates the internal orbit parameters.
    # @param xPos x-position
    # @param yPos y-position
    # @param xVel x-velocity
    # @param yVel y-velocity
    # @return e eccentricity vector
    # @return h specific orbital angular momentum
    ##
    def fromCartesianState (xPos, yPos, xVel, yVel):

        # Expand dimensionality of position and velocity (required for cross product)
        pos = np.array([xPos, yPos, 0])
        vel = np.array([xVel, yVel, 0])

        # Precompute the distance from the origin
        pos_norm = np.linalg.norm(pos)

        # Compute and update the orbital elements
        self.h = np.cross(pos, vel)
        self.e = np.cross(vel, self.h) / self.mu - pos / pos_norm
        return self.e, self.h


    ##
    # @brief TODO Computes kepler elementes from polar state vector.
    #        Updates the internal orbit parameters
    # @param rho
    # @param theta
    # @param rhoDot
    # @param thetaDot
    # @return e eccentricity vector
    # @return h specific orbital angular momentum
    ##
    def fromPolarState (rho, theta, rhoDot, thetaDot):
        # TODO
        return self.e, self.h


    ##
    # @brief Discretizes the orbit into a bunch of positions.
    #        Mainly used for plotting. Could also be used for
    #        pathfollowing (?). COULD BE OPTIMIZED.
    # @param N Number of samples
    # @return An Nx2 vector representing the discretized orbit
    ##
    def discretize (N):

        # Prepare samples container
        samples = np.zeros((N,2))

        # Precompute norms
        h_norm = np.linalg.norm(self.h)
        e_norm = np.linalg.norm(self.e)

        # Compute semi-latus rectum
        p = h_norm**2 / mu

        # Compute big semi-axis (a) and small semi-axis (b)
        a = p / (1 - e_norm**2)
        b = np.sqrt(p*a)

        # Compute ellipse midpoint coordinates
        pos_pe =   self.e * p / (e_norm * (1 + e_norm))
        pos_ap = - self.e * p / (e_norm * (1 - e_norm))
        pos_mid = (pos_ap + pos_pe) / 2.0

        # Compute rotation of the ellipse (and precompute sin, cos)
        beta = np.arctan2(self.e[1], self.e[0])
        sb = np.sin(beta)
        cb = np.cos(beta)

        # Discretize the ellipse
        angles = np.linspace(0, 2*np.pi, N)
        for k in range(N):
            sa = np.sin(angles[k])
            ca = np.cos(angles[k])

            samples[k,0] = pos_mid[0] + (a * ca * cb - b * sa * sb)
            samples[k,1] = pos_mid[1] + (a * ca * sb + b * sa * cb)

        # Return the discretized orbit
        return samples


## 
# Run this script to test the orbit functions
##
if __name__ == '__main__':
    print("TODO")