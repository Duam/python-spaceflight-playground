import numpy as np
from spaceflight_playground.constants import Universe
from spaceflight_playground.aux_models.kepler_orbit import KeplerOrbit

if __name__ == '__main__':
    # Import plotting library
    import matplotlib.pyplot as plt

    # Create an orbit instance
    orbit = KeplerOrbit()
    R = Universe.moon_radius

    # Print parameters
    print("Orbit parameters:")
    print("Eccentricity vector: " + str(orbit.eccentricity))
    print("Angular momentum vector: " + str(orbit.specific_angular_momentum))

    # Set the orbit using the orbital parameters directly
    e = np.array([0.3, 0.3, 0])
    h = np.array([0, 0, 3000000])

    orbit.setOrbitalElements(e, h)
    print(orbit)
    samples_0 = orbit.discretize()

    # Set the orbit using polar state vector
    orbit.fromPolarState(R + 20000.0, 0.0, 0.0, 0.00095046751314)
    print(orbit)
    samples_1 = orbit.discretize()

    # Set the orbit using ellipse parameters
    orbit.fromEllipseParams(0.1, 0.0, 1.0)
    print(orbit)
    samples_2 = orbit.discretize()

    orbit.fromCartesianState(0.0, R + 20000.0, 1500.0, 0.0)
    print(orbit)
    samples_3 = orbit.discretize()

    # Plot orbit
    # plt.plot(samples_0[:,0], samples_0[:,1])
    plt.plot(samples_1[:, 0], samples_1[:, 1])
    # plt.plot(samples_2[:,0], samples_2[:,1])
    plt.plot(samples_3[:, 0], samples_3[:, 1])

    # plt.xlim([-2,2])
    # plt.ylim([-2,2])
    plt.grid(True)
    plt.show()
