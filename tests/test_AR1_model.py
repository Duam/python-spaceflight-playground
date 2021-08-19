import numpy as np
import matplotlib.pyplot as plt
from spaceflight_playground.models.AR1_model import AR1_model

if __name__ == '__main__':
    ar1 = AR1_model(autoregressivity=0.3, mean=5, variance=10)
    print("AR1-model parameters:")
    print("Process mean = " + str(ar1.mean))
    print("Process variance = " + str(ar1.variance))
    print("Autoregression parameter phi = " + str(ar1.autoregressivity))
    print("Autoregression constant parameter c = " + str(ar1.c))
    print("Normal distribution stddev = " + str(ar1.normdist_stddev))
    print("Initialized sample value = " + str(ar1.xCur))

    N = 1000
    kAxis = range(N)
    samples = np.array([])
    for k in kAxis:
        samples = np.append(samples, ar1.update())

    plt.plot(kAxis, samples)
    plt.xlabel('Timestep')
    plt.ylabel('Sample value')
    plt.show()
