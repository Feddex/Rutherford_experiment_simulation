import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare

###this code is to verify if the random generation of initial positions is actually gaussian

# Use Agg backend for Matplotlib
import matplotlib
matplotlib.use('Agg')

# Define the new fit function
def fit_gaussian_first(x, xo, s):
    return (1 / np.sqrt(2*np.pi) * s) * np.exp(-((x-xo)**2) / (2 * s**2))

# Load data from file
data = np.loadtxt('dist_initial_pos.txt')

# Define histogram parameters
bins = 1200

# Create histogram
hist, bin_edges = np.histogram(data, bins=bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Initial parameter values for fitting
initial_params = [0., 1.44E-10]  # Example initial values for n and s

# Fit the histogram data
popt, pcov = curve_fit(fit_gaussian_first, bin_centers, hist, p0=initial_params)
fitted_xo, fitted_s = popt

# Calculate expected values and normalize them
expected_values = fit_gaussian_first(bin_centers, fitted_xo, fitted_s)
expected_values *= sum(hist) / sum(expected_values)  # Normalize to match the sum of observed values

# Calculate reduced chi-squared
chi2, p = chisquare(hist, f_exp=expected_values)
reduced_chi2 = chi2 / (len(hist) - 2)  # Adjust degrees of freedom

# Plot histogram and fit
plt.figure(figsize=(10, 6))
plt.hist(data, bins=bins, label='Data', alpha=0.75)
plt.plot(bin_centers, expected_values, 'r-', label='Fit: n=%.2e, s=%.2e' % (fitted_xo, fitted_s), linewidth=1)

# Add reduced chi-squared to the plot
plt.legend(loc='best')
plt.xlabel('Angle of deflection (radians)')
plt.ylabel('Occurrences')
plt.title('Histogram and Fit of Deflection Angles')
plt.text(1.5, max(hist) * 0.8, 'Reduced Chi-squared: %.2f' % reduced_chi2, fontsize=12, color='red')
plt.grid(True)

# Set y-axis limits
plt.ylim(0, max(hist))  # Adjust the factor 0.2 as needed to show less of the y-axis

# Show plot
plt.show()
