import numpy as np
import matplotlib.pyplot as plt

np.random.seed(9)

# Create variables
N = 100000  # Number of objects
t_clusters = [10]
t_born = []
t_ms = []
t_cool = []
dt = 1
alpha = 2.35
min_mass = 0.1
max_mass = 100
mass_ms = []
num_wd = 0
num_bh = 0
lum = []
t_eff = []
R = []
sigma = 5.67e-8  # W/m^2 K^4

# Monte Carlo Simulation
while len(mass_ms) < N:
    # Generate mass directly from the specified range
    rand_mass = np.random.uniform(min_mass, max_mass)

    # For mass distribution, just accept it directly
    mass_ms.append(rand_mass)

# Create T_ms and T_i array
for i in range(N):
    # Calculate t_ms
    t_msi = 10 * mass_ms[i] ** -3.5
    t_ms.append(t_msi)

    # Calculate t_born
    t_i = dt * np.random.uniform(0, 1)
    t_born.append(t_i)

# Finding t and luminosity etc
for t_cluster in t_clusters:
    for i in range(N):
        t_cool_i = t_cluster - t_born[i] - t_ms[i]
        t_cool.append(t_cool_i)

        # Calculate Luminosity
        if t_cool_i < 0:
            if mass_ms[i] >= 20:
                lum_i = 81 * mass_ms[i] ** 2.14
                lum.append(lum_i)
            elif 2 < mass_ms[i] <= 20:
                lum_i = 1.78 * mass_ms[i] ** 3.5
                lum.append(lum_i)
            elif mass_ms[i] <= 2:
                lum_i = 0.75 * mass_ms[i] ** 4.8
                lum.append(lum_i)

            # Calculate Radius R
            if mass_ms[i] >= 1.2:
                R_i = 10 ** 0.66 * np.log(mass_ms[i]) + 0.5
                R.append(R_i)
            elif mass_ms[i] < 1.2:
                R.append(mass_ms[i])

            # Calculate Temperature T_eff
            # Ensure that temperature calculation is correct
            t_eff_i = (lum[-1] / (4 * np.pi * R[-1] ** 2 * sigma)) ** (1 / 4)
            t_eff.append(t_eff_i)

        elif t_cool_i > 0:
            if mass_ms[i] < 10:
                num_wd += 1
            elif mass_ms[i] > 10:
                num_bh += 1

# Output the number of objects
print(f"Luminosity: {lum}")
print(f"Radiius: {R}")
print(f"Number of Black Holes or Neutron Stars: {num_bh}")
print(f"Number of White Dwarfs: {num_wd}")

# Plot results
plt.hist(mass_ms, bins=100, log=True, edgecolor='black')
plt.xlabel('Mass (Msun)')
plt.ylabel('Number of Stars')
plt.title('Mass Distribution')
plt.show()

plt.hist(t_born, bins=100, log=True, edgecolor='black')
plt.xlabel('t_born')
plt.ylabel('Number of Stars')
plt.title('TBORN')
plt.show()

# Scatter plot of Temperature vs Luminosity
plt.scatter(t_eff, lum, alpha=0.5)
plt.yscale('log')  # Set y-axis to log scale
# plt.xscale('log')  # Optionally, you can set x-axis to log scale too
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (W)')
plt.title('Plot of Temperature vs Luminosity')
plt.xlim(1, 2000)  # Limit x-axis for better visibility
plt.ylim(1e-5, 1e5)  # Set limits to match your desired scale
plt.show()
