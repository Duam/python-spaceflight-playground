#!/usr/bin/python3

import numpy as np

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

    ##
    # @brief Wrapper for numpy normal sampler
    ##
    def sample_normal(self):
        return np.random.normal(0, self.normdist_stddev)

    ## 
    # @brief Computes the next sample and updates the internal value
    # @return The next sample
    ##
    def update(self):
        self.xCur = self.c + self.phi * self.xCur + self.sample_normal()
        return self.xCur

    ## 
    # @brief Convenience method: Sample value getter
    # @return Current sample value
    ##
    def getCurrentSample(self):
        return self.xCur


##
# Execute this script to test the model
##    
if __name__ == '__main__':

    # Import plotting library
    import matplotlib.pyplot as plt

    # Create an AR1 model instance
    ar1 = AR1_model(phi=0.3, mean=5, variance=10)

    # Print out the parameters
    print("AR1-model parameters:")
    print("Process mean = " + str(ar1.mean))
    print("Process variance = " + str(ar1.variance))
    print("Autoregression parameter phi = " + str(ar1.phi))
    print("Autoregression constant parameter c = " + str(ar1.c))
    print("Normal distribution stddev = " + str(ar1.normdist_stddev))
    print("Initialized sample value = " + str(ar1.xCur))

    # Simulation parameters
    N = 1000
    kAxis = range(N)

    # Simulation result container
    samples = np.array([])

    # Simulate and store results
    for k in kAxis:
        samples = np.append(samples, ar1.update())

    # Plot results
    plt.plot(kAxis, samples)
    plt.xlabel('Timestep')
    plt.ylabel('Sample value')
    plt.show()