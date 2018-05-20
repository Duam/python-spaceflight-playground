# Model
This model is that of a two dimensional point mass rocket. Its thrusters can fire in any direction.

## Inputs
- u<sub>r<\sub>: Scaled thrust in radial direction (Between 0 and 1)
- u<sub>t<\sub>: Scaled thrust in tangential direction (Between 0 and 1)

## Disturbances
No disturbances

## States
- r: Altitude (km)
- &theta;: Angle from horizontal axis (microRad)
- rDot: Radial velocity (km per second)
- &theta;Dot: Angular velocity (microRad per second)
- m: Mass (kg)

# Assumptions/Simplifications:
- No atmosphere
- The moon is a point mass and perfectly circular
