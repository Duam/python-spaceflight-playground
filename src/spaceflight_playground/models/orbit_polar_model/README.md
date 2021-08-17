# Model
This model is that of a two dimensional point mass rocket. Its thrusters can fire in any direction.

## Inputs
- T_r: Scaled thrust in radial direction (Between 0 and 1)
- T_theta: Scaled thrust in tangential direction (Between 0 and 1)

## Disturbances
No disturbances

## States
- r: Altitude (km)
- theta: Angle from horizontal axis (microRad)
- rDot: Radial velocity (km per second)
- thetaDot: Angular velocity (microRad per second)
- m: Mass (kg)

# Assumptions/Simplifications:
- No atmosphere, because it's the moon
- The moon doesn't rotate
- Spacecraft has no rotation
- Thrusters can fire in any direction
- The moon is a point mass and perfectly circular
