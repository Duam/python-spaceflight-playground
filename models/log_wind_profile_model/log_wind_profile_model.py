#!/usr/bin/python3

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
    def __init__(self, zr=0.001, z0=10, u0=1):
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
        return u0 * np.log(z/zr) / np.log(z0/zr)