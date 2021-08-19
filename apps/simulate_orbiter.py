import numpy as np
import casadi as cas
import matplotlib.pyplot as plt
from spaceflight_playground.rk4step import rk4step_ode
from spaceflight_playground.aux_models.kepler_orbit import KeplerOrbit
from spaceflight_playground.orbiter.orbiter_model import (
    PolarOrbiter,
    PolarOrbiterState,
    PolarOrbiterThrust
)

spacecraft = PolarOrbiter()

# Simulation parameters
total_time = 600.0  # [s]
num_samples = 100
timestep = total_time / num_samples

# Create discretized system model with CasADi
x = cas.MX.sym('x', spacecraft.num_states, 1)
u = cas.MX.sym('u', spacecraft.num_forces, 1)
state_derivative_fun = cas.Function(
    'ode', [x, u],
    [spacecraft.ode_scaled(PolarOrbiterState(x), PolarOrbiterThrust(u)).vector],
    ['state', 'thrust'], ['state_derivative'])
next_state = rk4step_ode(state_derivative_fun, x, u, timestep)
next_state_fun = cas.Function('next_state', [x, u], [next_state], ['x', 'u'], ['xnext'])

# Initial state
initial_state = PolarOrbiterState(
    cas.DM([
        spacecraft.moon_radius,
        0,
        0,
        0,
        spacecraft.full_mass
    ])
).scale(spacecraft.scale)
print(f"Initial state: {initial_state.vector}")

# Choose controls for simulation
thrusts = np.zeros((spacecraft.num_forces, num_samples))
thrust_radial_stop = 20
thrust_angular_stop = 50
thrusts[0, 0:thrust_radial_stop] = 0.19 * np.ones(thrust_radial_stop)
thrusts[1, 0:thrust_angular_stop] = 0.36 * np.ones(thrust_angular_stop)

# Simulate the system
xs = cas.DM.zeros((spacecraft.num_states, num_samples + 1))
xs[:, 0] = initial_state.vector
for k in range(1, num_samples + 1):
    xs[:, k] = next_state_fun(xs[:, k - 1], thrusts[:, k - 1])

# Unscale trajectory
for k in range(num_samples + 1):
    xs[:, k] = xs[:, k] * spacecraft.unscale

xs = xs.full()

# Prepare plotting
tAxis = np.linspace(0, total_time - timestep, num_samples + 1)
plt.figure(1)

# Reference altitude [km] & angular velocity [µrad]
target_altitude = 20  # [km]
target_distance = spacecraft.moon_radius + target_altitude * 1e3
target_angular_velocity = np.sqrt(spacecraft.gravitational_constant * spacecraft.moon_mass / target_distance ** 3)  # [rad/s]

print("Reference altitude: " + str(target_altitude))
print("Reference angular velocity: " + str(target_angular_velocity))

# Plot
plt.subplot(321)
plt.plot(tAxis, xs[0, :])
plt.plot(tAxis, target_distance * np.ones((tAxis.size, 1)))
plt.ylabel('Altitude [km]')

plt.subplot(322)
plt.plot(tAxis, xs[1, :])
plt.plot(tAxis, initial_state.vector[1].full() * np.ones((tAxis.size, 1)))
plt.ylabel('Angle [murad]')

plt.subplot(323)
plt.plot(tAxis, xs[2, :])
plt.plot(tAxis, initial_state.vector[2].full() * np.ones((tAxis.size, 1)))
plt.ylabel('Radial vel. [km/s]')

plt.subplot(324)
plt.plot(tAxis, xs[3, :])
plt.plot(tAxis, target_angular_velocity * np.ones((tAxis.size, 1)))
plt.ylabel('Angular vel. [murad/s]')

plt.subplot(325)
plt.plot(tAxis, xs[4, :])
plt.ylabel('Mass [kg]')

plt.show()

quit(0)

# Convert trajectory to cartesian coordinates
for k in range(num_samples + 1):
    xs[0, k] = xs[0, k] + spacecraft.moon_radius

# Unscale trajectory
# for k in range(N+1):
#    xs[:,k] = xs[:,k] * spacecraft.unscale

xs_cart, us_cart = traj_pol2cart(xs[0:4, :], thrusts)

print(xs_cart[0, :])
print(xs_cart[1, :])

# Create an orbit instance
orb = KeplerOrbit()
orb.fromPolarState(
    distance=spacecraft.moon_radius + spacecraft.unscale[0] * target_altitude,
    angle=0.0,
    radial_velocity=0.0,
    angular_velocity=spacecraft.unscale[3] * target_angular_velocity
)

# Write trajectory to xml file
write_to_xml(
    filename='orbit_animator_trajectory.xml',
    T=total_time,
    N=num_samples,
    e=orb.eccentricity[0:2],
    h=orb.specific_angular_momentum[2],
    xPoses=xs_cart[0, :],
    yPoses=xs_cart[1, :],
    xVelos=xs_cart[2, :],
    yVelos=xs_cart[3, :],
    masses=xs[4, :],
    xForces=us_cart[0, :],
    yForces=us_cart[1, :]
)