# Dependencies
python3 with: casadi, numpy, matplotlib

# ToDo
* Orbit animator
* polar orbit model test with animation
* Check liftoff model equations (particularly torque and thrust)
* Spacecraft animator
* Passing around orbits and trajectories feels VERY uneasy and complicated. Maybe implement a trajectory class? Also extend xml interface that e.g. takes a kepler_orbit() instance as input and also outputs a kepler_orbit() instead of the orbital elements. Make xml writer its own class?
* Outsource unit tests to their own test folders?