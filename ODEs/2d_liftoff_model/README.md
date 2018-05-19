# Model
This model is that of a two dimensional rocket. It has one gimballed thruster at its
base.

There are three important points in space:
- Center of mass (COM): The gravitational force applies here
- Center of pressure (COP): The wind force applies here
- Base: Gimballed thruster force

## Inputs
- $F_T$ the thruster force
- $\mu$ the gimbal angle

## Disturbances
- $F_w$ the wind force

## States
- $\ţheta$ Angle, measured from the horizontal (inertial) plane
- $\dot{\theta}$ Angular velocity
- $x$, $y$ Position
- $\dot{x}$ $\dot{y}$ Linear velocit

# Assumptions:
- Wind direction stays horizontal
- Logarithmic wind profile
- Gravity points down
- 
