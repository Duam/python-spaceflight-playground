# Model
Logarithmic wind profile model. It models the winds magnitude along the vertical (height) axis.

## Initialization parameters
- zr: Surface roughness parameter (default=0.001)
- z0: Reference height (default=10)
- u0: Reference wind speed (default=1)

## Usage
1. Initialize the model using the above parameters
2. Use the getWindspeed(z) function to get the wind speed at an altitude of height z
