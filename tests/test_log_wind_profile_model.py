import numpy as np
import matplotlib.pyplot as plt
from spaceflight_playground.models.log_wind_profile_model import log_wind_profile_model

if __name__ == '__main__':
    # Create wind profile
    wind_profile = log_wind_profile_model()

    # Print parameters
    print("Wind profile parameters:")
    print("Surface roughness: " + str(wind_profile.zr) + " m")
    print("Reference height: " + str(wind_profile.z0) + " m")
    print("Reference wind speed: " + str(wind_profile.u0) + " m/s")

    # Create altitude axis
    zAxis = np.linspace(wind_profile.zr, 100.0, 100)

    # Compute wind profile along altitude
    windspeeds = np.zeros(zAxis.size)
    for k in range(zAxis.size):
        windspeeds[k] = wind_profile.getWindspeed(zAxis[k])

    # Plot wind profile
    plt.plot(windspeeds, zAxis)
    plt.ylabel("Altitude [m]")
    plt.xlabel("Wind speed [m/s]")
    plt.show()