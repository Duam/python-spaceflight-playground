# Model
This model is that of a two dimensional rocket. It has one gimballed thruster at its
base.

There are three important points in space:
- Center of mass (COM): The gravitational force applies here
- Center of pressure (COP): The wind force applies here
- Base: Gimballed thruster force

## Inputs
- F<sub>T<\sub>: the thruster force
- &mu;: the gimbal angle

## Disturbances
- F<sub>w<\sub>: the wind force

## States
- x, y: Position
- xdot, ydot Linear velocit
- &theta; Angle, measured from the vertical (inertial) plane
- &theta;dot$ Angular velocity

# Assumptions/Simplifications:
- Wind direction stays horizontal
- Logarithmic wind profile
- Gravity points down
