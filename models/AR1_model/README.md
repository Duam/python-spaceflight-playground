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