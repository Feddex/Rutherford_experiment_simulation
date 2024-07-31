import numpy as np
from scipy.integrate import quad

# Define the cotangent function
def cot(x):
    return np.cos(x) / np.sin(x)

# Define the integrand
def integrand(phi, a, sigma):
    return (a / (2 * np.sqrt(2 * np.pi * sigma**2))) * np.exp(-a**2 * (cot(phi / 2)**2) / (2 * sigma**2)) * (1 / np.sin(phi / 2)**2)

# Set the parameters
a = 1  # Example value for a
sigma = 1  # Example value for sigma

# Perform the integration over the interval [0, 2π]
I, error = quad(integrand, 0, 2 * np.pi, args=(a, sigma))

print(f"I: {I}, error: {error}")
