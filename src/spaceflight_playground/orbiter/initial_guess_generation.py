import casadi as cas
from typing import List

def triangle(
        num_samples: float,
        max_velocity: float,
) -> List[float]:
    """

    :param total_time:
    :param num_samples:
    :param max_velocity:
    :return:
    """
    half_samples = int(num_samples / 2.)
    slope = max_velocity / half_samples
    first_interval = [k * slope for k in range(half_samples)]
    second_interval = [ - slope * k + max_velocity for k in range(num_samples - half_samples)]
    return first_interval + second_interval

