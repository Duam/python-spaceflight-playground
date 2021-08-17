#!/usr/bin/python3

##
# @file main_orbit_polar.py
# @author Paul Daum
# @brief Optimize spacecraft trajectory to reach a circular orbit at 
#        a specified altitude. Uses multiple shooting.
##

import casadi as cas
import numpy as np
from spaceflight_playground.models.polar_orbiter import PolarOrbiter, PolarOrbiterState, PolarOrbiterThrust
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
integrator_substeps = 10
integrator_stepsize = timestep / integrator_substeps

# Values for terminal constraints
target_altitude = 20  # [km]
target_angular_velocity = 10 ** 3 * np.sqrt(spacecraft.gravitational_constant * spacecraft.moon_mass /\
                                            (spacecraft.moon_radius + 10 ** 3 * target_altitude) ** 3)

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

print(ode)

# Create an integrator for the ode
Xk = x
for k in range(integrator_substeps):
    Xk = rk4step_ode(ode, Xk, u, integrator_stepsize)
F = cas.Function('F', [x, u], [Xk], ['x','u'], ['xk'])

# Create stage cost for the OCP
l = u[0]**2 + u[1]**2
l = cas.Function('l', [x, u], [l], ['x','u'], ['l'])

# Create an integrator for the stage cost
Lk = 0
for k in range(integrator_substeps):
    Lk = rk4step_L(l, Lk, x, u, integrator_stepsize)
L = cas.Function('L', [x,u], [Lk], ['x','u'], ['L'])

# Create an initial guess for the OCP by forward simulation
us_init = np.zeros((num_samples, spacecraft.num_forces))
n_r_stop = 60
n_theta_stop = 85

us_r = np.zeros(num_samples)
us_r[0:n_r_stop] = 0.1075 * np.ones(n_r_stop)
us_init[:,0] = us_r

us_theta = np.zeros(num_samples)
us_theta[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
us_init[:,1] = us_theta

xs = cas.DM.zeros((num_samples + 1, spacecraft.num_states))
xs[0,:] = spacecraft.x0_scaled
for k in range(num_samples):
    xs[k+1,:] = F(xs[k,:],us_init[k,:])

xs_init = xs.full()

# Print debug message
print("== Initial guess ==")
print("xs_init size: " + str(xs_init.shape))
print("us_init size: " + str(us_init.shape))
print(xs_init)
print("Initial guess computed. Now starting creation of OCP.")

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
Xk = xs_init[0,:]
for k in range(num_samples):
    
    # NLP variable for control
    Uk = cas.MX.sym('U_' + str(k), spacecraft.num_forces, 1)
    w = cas.vertcat(w, Uk)
    lbw = cas.vertcat(lbw, -cas.inf, -cas.inf)
    ubw = cas.vertcat(ubw,  cas.inf,  cas.inf)
    w0 = cas.vertcat(w0, us_init[k,:])

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
        0.0,
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
    w0 = cas.vertcat(w0, xs_init[k+1,:])

    # Equality constraints to match intervals
    g = cas.vertcat(g, Xk_end - Xk)
    lbg = cas.vertcat(lbg, cas.DM.zeros(spacecraft.num_states, 1))
    ubg = cas.vertcat(ubg, cas.DM.zeros(spacecraft.num_states, 1))


# Terminal constraint on altitude
g = cas.vertcat(g, Xk_end[0])
lbg = cas.vertcat(lbg, target_altitude)
ubg = cas.vertcat(ubg, target_altitude)

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
print("w size: " + str(w.shape) + ", type: " + str(type(w)))
print("w0 size: " + str(w0.shape)  + ", type: " + str(type(w0)))
print("lbw size: " + str(lbw.shape) + ", type: " + str(type(lbw)))
print("ubw size: " + str(ubw.shape) + ", type: " + str(type(ubw)))
print("J size: " + str(J.shape) + ", type: " + str(type(J)))
print("g size: " + str(g.shape) + ", type: " + str(type(g)))
print("lbg size: " + str(lbg.shape)  + ", type: " + str(type(lbg)))
print("ubg size: " + str(ubg.shape) + ", type: " + str(type(ubg)))
print("Setting up and starting solver")

# Create an NLP solver
nlp = {}
nlp['f'] = J
nlp['x'] = w
nlp['g'] = g

opts = {}
#opts['ipopt.print_level'] = 0
opts['ipopt.print_info_string'] = 'yes'
solver = cas.nlpsol('solver', 'ipopt', nlp, opts)

# Solve the NLP
solver_in = {}
solver_in['x0'] = w0
solver_in['lbx'] = lbw
solver_in['ubx'] = ubw
solver_in['lbg'] = lbg
solver_in['ubg'] = ubg
solver_out = solver(**solver_in)
print("== OCP solved ==")

# Extract results
sol = solver_out['x']
print("sol size: " + str(sol.shape) + ", type: " + str(type(sol)))

u_opt = cas.DM.zeros((num_samples, spacecraft.num_forces))
x_opt = cas.DM.zeros((num_samples, spacecraft.num_states))

nxnu = spacecraft.num_states + spacecraft.num_forces

u_opt[:,0] = sol[0::nxnu]
u_opt[:,1] = sol[1::nxnu]
x_opt[:,0] = sol[2::nxnu]
x_opt[:,1] = sol[3::nxnu]
x_opt[:,2] = sol[4::nxnu]
x_opt[:,3] = sol[5::nxnu]
x_opt[:,4] = sol[6::nxnu]

# Write to .xml file