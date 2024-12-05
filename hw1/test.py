# This is some code for CPC that I was using to plot the Trinity Bomb shock wave.


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Data (convert time from ms to seconds)
t_data = np.array([0.10, 0.24, 0.38, 0.52, 0.66, 0.80, 0.94, 1.08, 1.22, 1.36, 1.50, 1.65, 1.79, 1.93, 3.26, 3.53, 3.80, 4.07, 4.34, 4.61, 15.0, 25.0, 34.0, 53.0, 62.0]) * 1e-3  # Convert ms to seconds
R_data = np.array([11.1, 19.9, 25.4, 28.8, 31.9, 34.2, 36.3, 38.9, 41.0, 42.8, 44.4, 46.0, 46.9, 48.7, 59.0, 61.1, 62.9, 64.3, 65.6, 67.3, 106.5, 130.0, 145.0, 175.0, 185.0])

# Logarithms of radius and time
log_R = np.log(R_data)
log_t = np.log(t_data)

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(log_t, log_R)

# Plot the data and the linear fit
plt.scatter(log_t, log_R, label='Data')
plt.plot(log_t, slope * log_t + intercept, color='red', label=f'Fit: slope={slope:.2f}')
plt.xlabel('log(t) [s]')
plt.ylabel('log(R) [m]')
plt.legend()
plt.show()

# Output slope to verify it corresponds to 2/5
print(f"Slope from linear regression: {slope:.2f}")

# Now calculate E (use slope from regression, which should approximate 2/5)
# Assume slope = 2/5, and calculate E
# Scaling law: R(t) ∝ E^(1/5) * t^(2/5)
# Calculate yield E based on the known proportionality
E = (1.28 * 1e3)**(5/3) * 36.2**5  # Simplified energy calculation
print(f"Estimated yield E (in Joules): {E:.2e}")
