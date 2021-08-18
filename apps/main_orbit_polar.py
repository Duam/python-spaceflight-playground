#!/usr/bin/python3

"""
:author: Paul Daum
:brief: Optimizes a spacecraft's trajectory to reach a circular orbit at a specified altitude.
    Uses multiple shooting.
"""

import casadi as cas
import numpy as np
import matplotlib.pyplot as plt
from spaceflight_playground.orbiter.orbiter_model import PolarOrbiter, PolarOrbiterState, PolarOrbiterThrust
from spaceflight_playground.rk4step import rk4step_L, rk4step_ode

# Create a spacecraft instance
spacecraft = PolarOrbiter()

# Print parameters
print("== Universe parameters ==")
print("Grav. const. G: " + str(spacecraft.gravitational_constant) + " m^3/(kg*s^2)")
print("Moon mass M: " + str(spacecraft.moon_mass) + " kg")
print("Moon radius R: " + str(spacecraft.moon_radius) + " m")
print("== Spacecraft parameters ==")
print("Initial mass: " + str(spacecraft.full_mass) + " kg")
print("Empty mass: " + str(spacecraft.empty_mass) + " kg")
print("Max. thrust: " + str(spacecraft.max_thrust) + " N")
print("Fuel consumption coeff:" + str(spacecraft.fuel_consumption_per_second) + " kg/s")
print("State scale vector: " + str(spacecraft.scale))
print("State unscale vector: " + str(spacecraft.unscale))


# Simulation parameters
total_time = 600.0
num_samples = 100
timestep = total_time / num_samples

# Integration parameters
integrator_substeps = 1
integrator_stepsize = timestep / integrator_substeps

# Values for terminal constraints
target_altitude = 20  # [km]
target_distance = spacecraft.moon_radius + 1e3 * target_altitude
target_angular_velocity = 1e3 * np.sqrt(spacecraft.gravitational_constant * spacecraft.moon_mass /
                                        target_distance ** 3)

# Print out values
print("== Terminal values ==")
print("Altitude: " + str(target_altitude) + " km")
print("Angular velocity: " + str(target_angular_velocity) + " microrad")

# Create system model in casadi
x = cas.MX.sym('x', spacecraft.num_states, 1)
u = cas.MX.sym('u', spacecraft.num_forces, 1)
state = PolarOrbiterState(x)
thrust = PolarOrbiterThrust(u)
state_derivative = spacecraft.ode_scaled(state, thrust).vector
ode = cas.Function('ode', [x, u], [spacecraft.ode_scaled(PolarOrbiterState(x), PolarOrbiterThrust(u)).vector],
                   ['state', 'thrust'], ['state_derivative'])

initial_state = PolarOrbiterState(
    cas.DM([
        spacecraft.moon_radius,
        0,
        0,
        0,
        spacecraft.full_mass
    ])).scale(spacecraft.scale)

print(f"Initial state: {initial_state.vector}")

# Create an integrator for the ode
Xk = x
for k in range(integrator_substeps):
    Xk = rk4step_ode(ode, Xk, u, integrator_stepsize)
F = cas.Function('F', [x, u], [Xk], ['x', 'u'], ['xk'])

# Create stage cost for the OCP
state_cost = u[0] ** 2 + u[1] ** 2
state_cost = cas.Function('l', [x, u], [state_cost], ['x', 'u'], ['l'])

# Create an integrator for the stage cost
Lk = 0
for k in range(integrator_substeps):
    Lk = rk4step_L(state_cost, Lk, x, u, integrator_stepsize)
L = cas.Function('L', [x, u], [Lk], ['x', 'u'], ['L'])

# Create an initial guess for the OCP by forward simulation
thrusts_initial_guess = np.zeros((num_samples, spacecraft.num_forces))
n_r_stop = 50
n_theta_stop = 50

radial_thrusts = np.zeros(num_samples)
radial_thrusts[0:n_r_stop] = 0.1075 * np.ones(n_r_stop)
thrusts_initial_guess[:, 0] = radial_thrusts

angular_thrusts = np.zeros(num_samples)
angular_thrusts[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
thrusts_initial_guess[:, 1] = angular_thrusts

states_initial_guess = cas.DM.zeros((num_samples + 1, spacecraft.num_states))
states_initial_guess[0, :] = initial_state.vector
for k in range(num_samples):
    states_initial_guess[k + 1, :] = F(states_initial_guess[k, :], thrusts_initial_guess[k, :])

states_initial_guess = states_initial_guess.full()

# Print debug message
print("== Initial guess ==")
print("xs_init size: " + str(states_initial_guess.shape))
print("us_init size: " + str(thrusts_initial_guess.shape))
print(states_initial_guess)
print("Initial guess computed. Now starting creation of OCP.")

opti = cas.Opti()
x0 = opti.parameter(spacecraft.num_states)
X = opti.variable(spacecraft.num_states, num_samples)
U = opti.variable(spacecraft.num_forces, num_samples - 1)

opti.minimize(sum([L(X[:, k], U[:, k]) for k in range(num_samples - 1)]))
opti.subject_to(X[:, 0] == x0)
opti.subject_to([X[:, k+1] == F(X[:, k], U[:, k]) for k in range(num_samples - 1)])
opti.subject_to([U[0, k]**2 + U[1, k]**2 <= 1 for k in range(num_samples - 1)])
opti.subject_to(X[0, -1] == target_distance)
opti.subject_to(X[2, -1] == 0)
opti.subject_to(X[3, -1] == target_angular_velocity)

opti.set_value(x0, initial_state.vector)
opti.set_initial(X, states_initial_guess[1:, :].T)
opti.set_initial(U, thrusts_initial_guess[1:, :].T)
opti.solver('ipopt', {'expand': True})
solution = opti.solve()

Xsol = solution.value(X)


fig, axs = plt.subplots(5, 1)
plt.sca(axs[0])
plt.plot(Xsol[0, :])
plt.sca(axs[1])
plt.plot(Xsol[1, :])
plt.sca(axs[2])
plt.plot(Xsol[2, :])
plt.sca(axs[3])
plt.plot(Xsol[3, :])
plt.sca(axs[4])
plt.plot(Xsol[4, :])

plt.show()
print(f"")
print(solution.value(X))

quit(0)
# Create the optimization variables
X = cas.MX.sym('X', spacecraft.num_states, num_samples)
U = cas.MX.sym('U', spacecraft.num_forces, num_samples)

# Start with empty NLP
w = []      # Optimization variables (xs, us)
w0 = []     # Initial guess
lbw = []    # Lower bound on opt. variables
ubw = []    # Upper bound on opt. variables
J = 0       # Cost function
g = []      # Nonlinear constraints
lbg = []    # Lower bound on constraints
ubg = []    # Upper bound on constraints

# Formulate NLP
Xk = states_initial_guess[0, :]
for k in range(num_samples):
    
    # NLP variable for control
    Uk = cas.MX.sym('U_' + str(k), spacecraft.num_forces, 1)
    w = cas.vertcat(w, Uk)
    lbw = cas.vertcat(lbw, -cas.inf, -cas.inf)
    ubw = cas.vertcat(ubw,  cas.inf,  cas.inf)
    w0 = cas.vertcat(w0, thrusts_initial_guess[k, :])

    # Circle constraints on controls
    g = cas.vertcat(g, Uk[0]**2 + Uk[1]**2)
    lbg = cas.vertcat(lbg, 0)
    ubg = cas.vertcat(ubg, 1)

    # Integrate till the end of the interval
    Xk_end = F(Xk, Uk)
    J = J + L(Xk, Uk)

    # New NLP variable for state
    Xk = cas.MX.sym('X_' + str(k+1), spacecraft.num_states, 1)
    w = cas.vertcat(w, Xk)
    lbw_k = cas.vertcat(
        spacecraft.moon_radius * spacecraft.scale[0],
        -cas.inf,
        -cas.inf,
        -cas.inf,
        spacecraft.empty_mass * spacecraft.scale[4]
    )
    ubw_k = cas.vertcat(
        cas.inf,
        cas.inf,
        cas.inf,
        cas.inf,
        spacecraft.full_mass * spacecraft.scale[4]
    )
    lbw = cas.vertcat(lbw, lbw_k)
    ubw = cas.vertcat(ubw, ubw_k)
    w0 = cas.vertcat(w0, states_initial_guess[k+1, :])

    # Equality constraints to match intervals
    g = cas.vertcat(g, Xk_end - Xk)
    lbg = cas.vertcat(lbg, cas.DM.zeros(spacecraft.num_states, 1))
    ubg = cas.vertcat(ubg, cas.DM.zeros(spacecraft.num_states, 1))


# Terminal constraint on distance state
g = cas.vertcat(g, Xk_end[0])
lbg = cas.vertcat(lbg, target_distance)
ubg = cas.vertcat(ubg, target_distance)

# Terminal constraint on radial velocity
g = cas.vertcat(g, Xk_end[2])
lbg = cas.vertcat(lbg, 0)
ubg = cas.vertcat(ubg, 0)

# Terminal constraint on angular velocity
g = cas.vertcat(g, Xk_end[3])
lbg = cas.vertcat(lbg, target_angular_velocity)
ubg = cas.vertcat(ubg, target_angular_velocity)

# Print debug message
print("== OCP created ==")
print(f"w size: {w.shape}, type: {type(w)}")
print(f"w0 size: {w0.shape}, type: {type(w0)}")
print(f"lbw size: {lbw.shape}, type: {type(lbw)}")
print(f"ubw size: {ubw.shape}, type: {type(ubw)}")
print(f"J size: {J.shape}, type: {type(J)}")
print(f"g size: {g.shape}, type: {type(g)}")
print(f"lbg size: {lbg.shape}, type: {type(lbg)}")
print(f"ubg size: {ubg.shape}, type: {type(ubg)}")
print(f"Setting up and starting solver")

# Create an NLP solver and solve it
nlp = {'f': J, 'x': w, 'g': g}
opts = {'expand': True, 'ipopt.print_info_string': 'yes'}
solver = cas.nlpsol('solver', 'ipopt', nlp, opts)
solver_in = {'x0': w0, 'lbx': lbw, 'ubx': ubw, 'lbg': lbg, 'ubg': ubg}
solver_out = solver(**solver_in)
print("== OCP solved ==")

# Extract results
sol = solver_out['x']
print("sol size: " + str(sol.shape) + ", type: " + str(type(sol)))

u_opt = cas.DM.zeros((num_samples, spacecraft.num_forces))
x_opt = cas.DM.zeros((num_samples, spacecraft.num_states))

nxnu = spacecraft.num_states + spacecraft.num_forces

u_opt[:, 0] = sol[0::nxnu]
u_opt[:, 1] = sol[1::nxnu]
x_opt[:, 0] = sol[2::nxnu]
x_opt[:, 1] = sol[3::nxnu]
x_opt[:, 2] = sol[4::nxnu]
x_opt[:, 3] = sol[5::nxnu]
x_opt[:, 4] = sol[6::nxnu]

# Write to .xml file