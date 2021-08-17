#!/usr/bin/python3

## 
# @file log_wind_profile_model.py
# @author Paul Daum
##

import numpy as np

## 
# @class log_wind_profile_model
# @brief A logarithmic wind profile model.
##
class log_wind_profile_model:

    ##
    # @brief Initialization procedure
    # @param zr Surface roughness parameter.
    #           Can be found in tables online.
    # @param z0 Reference height
    # @param u0 Reference wind speed
    ##
    def __init__(self, zr=0.001, z0=10.0, u0=1.0):
        # Surface roughtness parameter
        self.zr = zr
        # Reference height
        self.z0 = z0
        # Reference wind speed
        self.u0 = u0

    ##
    # @brief Computes the the wind speed according 
    #        to the logarithmic wind profile model.
    # @param z The height of which we want to get
    #          the wind speed
    # @return The wind speed
    ##
    def getWindspeed(self, z):
        # Check for invalid input
        if (z == 0):
            print("In log_wind_profile_model: Input z can not be zero. Check your input.")
            return 0

        # Compute wind speed and return
        u = self.u0 * np.log(z/self.zr) / np.log(self.z0/self.zr)
        return u


##
# Run this script to test the model
##
if __name__ == '__main__':
    
    # Import plotting library
    import matplotlib.pyplot as plt
    
    # Create wind profile 
    wind_profile = log_wind_profile_model()

    # Print parameters
    print("Wind profile parameters:")
    print("Surface roughness: " + str(wind_profile.zr) + " m")
    print("Reference height: " + str(wind_profile.z0) + " m")
    print("Reference wind speed: " + str(wind_profile.u0) + " m/s")

    # Create altitude axis
    zAxis = np.linspace(wind_profile.zr, 100.0, 100)

    # Compute wind profile along altitude
    windspeeds = np.zeros(zAxis.size)
    for k in range(zAxis.size):
        windspeeds[k] = wind_profile.getWindspeed(zAxis[k])

    # Plot wind profile
    plt.plot(windspeeds, zAxis)
    plt.ylabel("Altitude [m]")
    plt.xlabel("Wind speed [m/s]")
    plt.show()