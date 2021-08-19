import numpy as np


def compute_log_wind_speed(
        height: float,
        surface_roughness: float,
        reference_height: float,
        reference_wind_speed: float
) -> float:
    """Computes the wind speed according to a logarithmic wind profile model. It models
    the wind's magnitude along the vertical (height) axis.
    :param height: The height to compute the wind speed at [m].
    :param surface_roughness: The location's surface roughness. Can be found in tables online.
    :param reference_height: The height at which reference_wind_speed is measured.
    :param reference_wind_speed: The wind speed at the reference_height.
    :return: The wind speed at the given height.
    """
    assert height >= 0, f"Height must be non-negative, is {height}m."
    return reference_wind_speed * np.log(height / surface_roughness) / np.log(reference_height / surface_roughness)
    u = self.u0 * np.log(z / self.zr) / np.log(self.z0 / self.zr)
