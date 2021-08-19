import numpy as np
import casadi as cas
import matplotlib.pyplot as plt
from spaceflight_playground.rk4step import rk4step_ode
from spaceflight_playground.aux_models.kepler_orbit import KeplerOrbit
from spaceflight_playground.orbiter.orbiter_model import PolarOrbiter
from spaceflight_playground.conversion import *
from spaceflight_playground.xml_writer import *

if __name__ == '__main__':
    spacecraft = PolarOrbiter()

    # Simulation parameters
    T = 600.0
    N = 100
    DT = T / N

    # Integration parameters
    nn = 10  # Integrator steps per step
    h = DT / nn  # Step size of integrator step

    # Create system model with CasADi
    x = cas.MX.sym('x', spacecraft.num_states, 1)
    u = cas.MX.sym('u', spacecraft.num_forces, 1)
    f = cas.Function('f', [x, u], [spacecraft.ode_scaled(x, u)])

    # Discretize spacecraft dynamics using rk4
    Xk = x
    for k in range(nn):
        Xk = rk4step_ode(f, Xk, u, h)

    F = cas.Function('F', [x, u], [Xk], ['x', 'u'], ['xk'])

    # Choose controls for simulation
    us = np.zeros((spacecraft.num_forces, N))
    n_r_stop = 60
    n_theta_stop = 85

    us_r = np.zeros(N)
    us_r[0:n_r_stop] = 0.1075 * np.ones(n_r_stop)
    us[0, :] = us_r

    us_theta = np.zeros(N)
    us_theta[0:n_theta_stop] = 0.2 * np.ones(n_theta_stop)
    us[1, :] = us_theta

    # Simulate the system
    xs = cas.DM.zeros((spacecraft.num_states, N + 1))
    xs[:, 0] = spacecraft.x0_scaled

    for k in range(1, N + 1):
        xs[:, k] = F(xs[:, k - 1], us[:, k - 1])

    # Unscale trajectory
    for k in range(N + 1):
        xs[:, k] = xs[:, k] * spacecraft.unscale

    xs = xs.full()

    # Prepare plotting
    tAxis = np.linspace(0, T - DT, N + 1)
    plt.figure(1)

    # Reference altitude [km] & angular velocity [µrad]
    altitude_ref = 20
    angVel_ref = 1e3 * np.sqrt(spacecraft.mu / (spacecraft.moon_radius + 1e3 * altitude_ref) ** 3)

    print("Reference altitude: " + str(altitude_ref))
    print("Reference angular velocity: " + str(angVel_ref))

    # Plot
    plt.subplot(321)
    plt.plot(tAxis, xs[0, :])
    plt.plot(tAxis, altitude_ref * np.ones((tAxis.size, 1)))
    plt.ylabel('Altitude [km]')

    plt.subplot(322)
    plt.plot(tAxis, xs[1, :])
    plt.plot(tAxis, spacecraft.x0[1].full() * np.ones((tAxis.size, 1)))
    plt.ylabel('Angle [murad]')

    plt.subplot(323)
    plt.plot(tAxis, xs[2, :])
    plt.plot(tAxis, spacecraft.x0[2].full() * np.ones((tAxis.size, 1)))
    plt.ylabel('Radial vel. [km/s]')

    plt.subplot(324)
    plt.plot(tAxis, xs[3, :])
    plt.plot(tAxis, angVel_ref * np.ones((tAxis.size, 1)))
    plt.ylabel('Angular vel. [murad/s]')

    plt.subplot(325)
    plt.plot(tAxis, xs[4, :])
    plt.ylabel('Mass [kg]')

    plt.show()

    # Convert trajectory to cartesian coordinates
    for k in range(N + 1):
        xs[0, k] = xs[0, k] + spacecraft.moon_radius

    # Unscale trajectory
    # for k in range(N+1):
    #    xs[:,k] = xs[:,k] * spacecraft.unscale

    xs_cart, us_cart = traj_pol2cart(xs[0:4, :], us)

    print(xs_cart[0, :])
    print(xs_cart[1, :])

    # Create an orbit instance
    orb = KeplerOrbit()
    orb.fromPolarState(
        distance=spacecraft.moon_radius + spacecraft.unscale[0] * altitude_ref,
        angle=0.0,
        radial_velocity=0.0,
        angular_velocity=spacecraft.unscale[3] * angVel_ref
    )

    # Write trajectory to xml file
    write_to_xml(
        filename='orbit_animator_trajectory.xml',
        T=T,
        N=N,
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