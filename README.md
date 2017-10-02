# ControlledRocket

Here I experiment with different rocket models that are controlled using optimal control techniques.
The models will be regularly expanded and different control techniques will be tested for that purpose.

## V0: Point mass rocket in 2D
This is a rocket which is modelled as a point mass, which lifts off of the surface of the 
moon, which is also modelled as a point mass. Its target is a stable, circular orbit with
a given radius. The rocket does not have a rotation and its thrusters can fire in any 
arbitrary direction.
### States
x(1) = r        -- Distance from the origin (km) <br \>
x(2) = theta    -- Angle from horizontal axis (milliRad) <br \>
x(3) = rDot     -- Radial velocity (km per second) <br \>
x(4) = thetaDot -- Angular velocity (milliRad per second) <br \>
x(5) = m        -- Mass (kg) <br \>
### Controls
u(1) = u_r      -- Thrust in radial direction (percentage of max. thrust) <br \>
u(2) = u_t      -- Thrust in tangential direction (percentage of max. thrust) <br \>

