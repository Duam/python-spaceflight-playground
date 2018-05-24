#!/usr/bin/python3

##
# @file kepler_orbit_2d.py
# @author Paul Daum
##

import numpy as np

import sys, os
sys.path.append(os.path.realpath('../../'))
sys.path.append(os.getcwd())

from utils.conversion import pol2cart, state_pol2cart

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

        # Eccentricity vector (in 2d only the first two elements are non-zero)
        self.e = np.array([0,0,0])
        # Specific orbital angular momentum (in 2d only the last element is non-zero)
        self.h = np.array([0,0,0])
        
    ##
    # @brief Returns a string that contains the orbital elements
    ##
    def toString(self):
        return "e=" + str(self.e) + ", h=" + str(self.h)


    ##
    # @brief Setter for the orbital elements
    # @param e Eccentricity vector (3x1)
    # @param h Specific orbital angular momentum vector (3x1)
    # @return Nothing
    ##
    def setOrbitalElements(self, e, h):
        self.e = e
        self.h = h

    ##
    # @brief Computes kepler elements from common ellipse parameters.
    # @param e Eccentricity (Scalar value)
    # @param angle Angle of the orbit
    # @param a Big semi-major axis
    # @return e Eccentricty vector
    # @return h Specific orbital angular momentum
    ## 
    def fromEllipseParams(self, e, angle, a):
        
        # Compute semi-latus rectup
        p = a * (1 - e**2)

        # Compute norm of specific orbital angular momentum
        h = np.sqrt(p * self.mu)

        # Since we are in 2d, that's enough for to momentum vector
        self.h = np.array([0, 0, h])

        # Compute elements of the eccentricity vector
        e_x = e * np.cos(angle)
        e_y = e * np.sin(angle)

        # Since we are in 2d, that's enought for the eccentricity vector
        self.e = np.array([e_x, e_y, 0])

        # return
        return self.e, self.h


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
    def fromCartesianState (self, xPos, yPos, xVel, yVel):

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
    # @brief Computes kepler elementes from polar state vector.
    #        Updates the internal orbit parameters. TODO COULD BE OPTIMIZED.
    # @param rho
    # @param theta
    # @param rhoDot
    # @param thetaDot
    # @return e eccentricity vector
    # @return h specific orbital angular momentum
    ##
    def fromPolarState (self, rho, theta, rhoDot, thetaDot):

        # Convert polar coordinates to cartesian coordinates
        x_pol = np.array([rho, theta, rhoDot, thetaDot])
        x_cart = state_pol2cart(x_pol)
        
        # Get the vector elements
        xPos = x_cart[0]
        yPos = x_cart[1]
        xVel = x_cart[2]
        yVel = x_cart[3]

        # Compute orbital elements and return
        self.fromCartesianState(xPos, yPos, xVel, yVel)
        return self.e, self.h


    ##
    # @brief Discretizes the orbit into a bunch of positions.
    #        Mainly used for plotting. Could also be used for
    #        pathfollowing (?).
    # @param N Number of samples
    # @return An Nx2 vector representing the discretized orbit
    ##
    def discretize (self, N=360):

        # Precompute norms
        e_norm = np.linalg.norm(self.e)
        h_norm = np.linalg.norm(self.h)

        # Compute semi-latus rectum
        p = h_norm**2 / self.mu

        # Compute rotation of the ellipse
        beta = np.arctan2(self.e[1], self.e[0])
        
        # Prepare samples container
        samples = np.zeros((N,2))

        # Discretize the ellipse
        angles = np.linspace(0, 2*np.pi, N)
        for k in range(N):
            # Get current angle
            theta = angles[k]
            
            # Compute radius at current angle ("orbit equation")
            rho = p / (1 + e_norm * np.cos(theta))

            # Rotate the orbit clockwise
            theta = theta + beta

            # Convert into cartesian coordinates
            pos = pol2cart(theta, rho)

            # Store coordinates in container
            samples[k,0] = pos[0]
            samples[k,1] = pos[1]

        # Return the discretized orbit
        return samples


## 
# Run this script to test the orbit functions
##
if __name__ == '__main__':

    # Import plotting library
    import matplotlib.pyplot as plt

    # Create an orbit instance
    orbit = kepler_orbit_2d()

    # Print parameters
    print("Orbit parameters:")
    print("Eccentricity vector: " + str(orbit.e))
    print("Angular momentum vector: " + str(orbit.h))

    # Set the orbit using the orbital parameters directly
    e = np.array([0.3, 0.3, 0])
    h = np.array([0, 0, 3000000])

    orbit.setOrbitalElements(e,h)
    print(orbit.toString())
    samples_0 = orbit.discretize()

    # Set the orbit using polar state vector
    orbit.fromPolarState(orbit.R + 20000.0, 0.0, 0.0, 0.00095046751314)
    print(orbit.toString())
    samples_1 = orbit.discretize()

    # Set the orbit using ellipse parameters
    orbit.fromEllipseParams(0.1, 0.0, 1.0)
    print(orbit.toString())
    samples_2 = orbit.discretize()
    
    # Plot orbit
    plt.plot(samples_0[:,0], samples_0[:,1])
    plt.plot(samples_1[:,0], samples_1[:,1])
    plt.plot(samples_2[:,0], samples_2[:,1])

    #plt.xlim([-2,2])
    #plt.ylim([-2,2])
    plt.grid(True)
    plt.show()