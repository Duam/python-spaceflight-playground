import numpy as np

"""
# Model
Logarithmic wind profile model. It models the winds magnitude along the vertical (height) axis.

## Initialization parameters
- zr: Surface roughness parameter (default=0.001)
- z0: Reference height (default=10)
- u0: Reference wind speed (default=1)

## Usage
1. Initialize the model using the above parameters
2. Use the getWindspeed(z) function to get the wind speed at an altitude of height z
"""

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
