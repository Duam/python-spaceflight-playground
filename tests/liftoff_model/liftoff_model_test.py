#!/usr/bin/python3

## 
# @file orbit_animator_test.py
# @author Paul Daum
##

import sys, os
sys.path.append(os.path.realpath('../../'))
sys.path.append(os.getcwd())

import numpy as np
import matplotlib.pyplot as plt
import casadi as cas

from integrators.rk4step import rk4step_ode
from models.liftoff_model.liftoff_model import liftoff_model
from models.liftoff_model.liftoff_trajectory import liftoff_trajectory


# Create a spacecraft instance
spacecraft = liftoff_model()

# Print some parameters
print("Spacecraft parameters:")
print("Grav. accel. g: " + str(spacecraft.g) + " m/s^2")
print("Distance base to COM: " + str(spacecraft.L) + " m")
print("Distance base to COP: " + str(spacecraft.l) + " m")
print("Mass: " + str(spacecraft.m) + " kg")
print("Max. thrust: " + str(spacecraft.maxThrust) + " N")
print("Initial state: " + str(spacecraft.x0))
print("Number of states: " + str(spacecraft.nx))
print("Number of disturbances: " + str(spacecraft.nd))
print("Number of controls: " + str(spacecraft.nu))

# Simulation parameters
T = 600.0
N = 100
DT = T/N

# Integration parameters 
nn = 10     # Integrator steps per step
h = DT/nn   # Step size of integrator step

# Create system model with CasADi
x = cas.MX.sym('x', spacecraft.nx, 1)
u = cas.MX.sym('u', spacecraft.nu + spacecraft.nd, 1)
f = cas.Function('f', [x,u], [spacecraft.dynamics(x,u)])

# Discretize spacecraft dynamics using rk4
Xk = x
for k in range(nn):
    Xk = rk4step_ode(f, Xk, u, h)

F = cas.Function('F', [x,u], [Xk], ['x','u'], ['xk'])

# Choose controls for simulation
wind_forces = np.zeros((spacecraft.nd, N))
controls = np.zeros((spacecraft.nu, N))
controls[0,:] = 0.7 * np.ones(N)
controls[1,50:51] = 0.000001 * np.ones(1)

us = np.append(controls, wind_forces, axis=0)

# Test evaluation of the ODE
print(us[:,0])
xdot_test = spacecraft.dynamics(spacecraft.x0, us[:,0])
print(xdot_test)

# Simulate the system
xs = cas.DM.zeros((spacecraft.nx, N+1))
xs[:,0] = spacecraft.x0

for k in range(1,N+1):
    xs[:,k] = F(xs[:,k-1], us[:,k-1])

xs = xs.full()

# Prepare plotting
tAxis = np.linspace(0, T, N+1)
plt.figure(1)

# Plot
plt.subplot(321)
plt.plot(tAxis, xs[0,:])
plt.ylabel('x-pos [m]')

plt.subplot(322)
plt.plot(tAxis, xs[1,:])
plt.ylabel('y-pos [m]')

plt.subplot(323)
plt.plot(tAxis, xs[2,:])
plt.ylabel('x-vel [m/s]')

plt.subplot(324)
plt.plot(tAxis, xs[3,:])
plt.ylabel('y-vel [m/s]')

plt.subplot(325)
plt.plot(tAxis, xs[4,:])
plt.ylabel('angle [rad]')

plt.subplot(326)
plt.plot(tAxis, xs[5,:])
plt.ylabel('angular vel. [rad/s]')

plt.show()


trajectory = liftoff_trajectory(T, N, spacecraft)
trajectory.setXs(xs)
trajectory.setUs(controls)
trajectory.setDs(wind_forces)

trajectory.toXML('trajectory.xml')