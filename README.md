# Dependencies
python3 with: casadi, numpy, matplotlib

# ToDo
1. Liftoff animator
2. Implement orbit trajectory class, change orbit_polar_model accordingly
3. Solve an OCP for liftoff

# Folder structure
|-- integrators
|-- ocps: Optimal control problems
|-- legacy: Old matlab code
|-- models
|   |-- AR1_model: Autoregression model for wind simulation
|   |-- kepler_orbit: 2D kepler orbit model
|   |-- liftoff_model: 2D spacecraft model at liftoff
|   |-- log_wind_profile_model: Logarithmic wind profile model
|   |-- orbit_cartesian_model: 2D spacecraft model in orbit (cartesian coordinates)
|   |-- orbit_polar_model: 2D spacecraft model in orbit (polar coordinates)
|-- tests
|   |-- liftoff_animator: Test for liftoff animator (TODO: merge with test for liftoff model)
|   |-- liftoff_model: Test for liftoff model
|   |-- orbit_animator: Test for orbit animator (TODO: merge with test for orbit model)
|-- utils