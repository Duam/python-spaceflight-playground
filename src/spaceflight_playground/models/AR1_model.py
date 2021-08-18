#!/usr/bin/python3

import numpy as np

"""
# Model
This model is used for simulating windspeeds in one direction. The wind
changes each second, therefore it's modelled using a gaussian distribution.
The timeseries is generated using an AR(1) model.

## Initialization parameters
- phi: The AR(1) parameter (default = 0.3)
- mean: Mean value of the timeseries (default = 0)
- variance: Variance of the timeseries (default = 1)

## Usage
1. Initialize the model using the above parameters
2. Use the update() function to compute the next sample. It also returns the new sample

- Use getCurrentSample() to get the current sample without triggering the update.
- Use updateParameters(phi, mean, variance) to update the parameters whenever they change
"""

##
# @class AR1_model
# @brief An order 1 autoregressive model with
#        a gaussian seed. Used for simulating
#        wind speeds in one direction.
##
class AR1_model:

    ##
    # @brief Initialization procedure
    # @param phi Autoregression parameter
    # @param mean Mean of the timeseries
    # @param variance Variance of the timeseries (must be > 0)
    ##
    def __init__(self, phi=0.3, mean=0, variance=1):
        # Process mean
        self.mean = mean
        # Process variance
        self.variance = variance
        # Autoregression parameter
        self.phi = phi
        # Constant term in the AR(1) process (derived parameter)
        self.c = self.mean * (1 - self.phi)
        # Standard deviation of the gaussian seed (derived parameter)
        self.normdist_stddev = np.sqrt((1 - self.phi**2) * self.variance)
        # Current sample initialization
        self.xCur = mean
        
    ##
    # @brief Parameter update function
    # @param phi Autoregression parameter
    # @param mean Mean of the timeseries
    # @param variance Variance of the timeseries (must be > 0)
    ##
    def updateParameters(self, phi=0.3, mean=0, variance=1):
        self.mean = mean
        self.variance = variance
        self.phi = phi
        self.c = self.mean * (1 - self.phi)
        self.normdist_stddev = np.sqrt((1 - self.phi**2) * self.variance)


    def sample_normal(self):
        """Wrapper for numpy normal sampler
        :return: Normally distributed sample.
        """
        return np.random.normal(0, self.normdist_stddev)


    def update(self):
        """Computes the next sample and updates the internal value.
        :return: The next sample.
        """
        self.xCur = self.c + self.phi * self.xCur + self.sample_normal()
        return self.xCur


    def getCurrentSample(self):
        """Convenience method: sample value getter.
        :return: Current sample value.
        """
        return self.xCur

