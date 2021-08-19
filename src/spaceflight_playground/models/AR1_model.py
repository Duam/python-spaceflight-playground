import numpy as np

class AR1_model:
    def __init__(self, autoregressivity: float=0.3, mean: float=0, variance: float=1):
        """An order 1 autoregressive model with a gaussian seed. Used for simulating
        wind speeds in one direction.
        :param autoregressivity: Autoregression parameter. Higher -> values change slower.
        :param mean: The mean of the timeseries.
        :param variance: The variance of timeseries.
        """
        self.set_parameters(autoregressivity, mean, variance)
        self.xCur = mean


    def set_parameters(self, autoregressivity: float, mean: float, variance: float) -> None:
        """Parameter update function.
        :param autoregressivity: Autoregression parameter. Higher -> values change slower.
        :param mean: Mean of the timeseries.
        :param variance: Variance of the timeseries.
        :return: None
        """
        assert variance > 0, f"Variance must be positive, ist {variance}"
        self.mean = mean
        self.variance = variance
        self.autoregressivity = autoregressivity


    def sample_normal(self):
        """Wrapper for numpy normal sampler
        :return: Normally distributed sample.
        """
        return np.random.normal(0, np.sqrt((1 - self.autoregressivity ** 2) * self.variance))


    def update(self):
        """Computes the next sample and updates the internal value.
        :return: The next sample.
        """
        self.xCur = self.mean * (1 - self.autoregressivity) + self.autoregressivity * self.xCur + self.sample_normal()
        return self.xCur


    def getCurrentSample(self):
        """Convenience method: sample value getter.
        :return: Current sample value.
        """
        return self.xCur

