import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(9)

# Constants
N = 1000           # Number of objects
cluster_ages = [8, 10, 12]  # Cluster ages in Gyrs
dt = 1             # Timestep
alpha = 2.35       # IMF constant
min_mass = 0.1     # Minimum star mass
max_mass = 100     # Maximum star mass
sigma = 5.67e-8    # W/m^2 K^4
C = 0.0101         # Constant
R_solar = 6.969e8  # Solar radius (meters)
L_solar = 3.827e26 # Solar luminosity (watts)

# To store results for all cluster ages
all_lum_ms = []  # Main sequence luminosities
all_t_eff = []   # Main sequence effective temperatures
all_lum_wd = []  # White dwarf luminosities
all_t_eff_wd = []  # White dwarf effective temperatures

# Loop over each cluster age
for t_cluster in cluster_ages:
    # Initialize arrays
    mass_ms, t_ms, t_born, t_cool = [], [], [], []
    lum_ms, lum_wd, t_eff, t_eff_wd = [], [], [], []
    num_wd, num_bh = 0, 0  # Counters for white dwarfs and black holes

    # Generate Main Sequence Masses using Monte Carlo
    while len(mass_ms) < N:
        rand_y = np.random.uniform(0, 1)
        rand_mass = np.random.uniform(min_mass, max_mass)
        mass_i = rand_mass**(-alpha)
        if rand_y < mass_i:
            mass_ms.append(rand_mass)

    # Calculate t_ms, t_born, and t_cool for each star
    for i in range(N):
        t_msi = 10 * mass_ms[i] ** -3.5
        t_ms.append(t_msi)
        t_i = dt * np.random.uniform(0, 1)
        t_born.append(t_i)
        t_cool_i = t_cluster - t_born[i] - t_ms[i]
        t_cool.append(t_cool_i)

    # Determine luminosity and temperature
    for i in range(N):
        mass = mass_ms[i]
        if t_cool[i] <= 0:  # Main Sequence
            if mass > 55:
                lum_i = 32000 * mass
            elif 2 < mass <= 55:
                lum_i = 1.4 * mass ** 3.5
            elif 0.43 < mass <= 2:
                lum_i = mass ** 4
            else:
                lum_i = 0.23 * mass ** 2.3
            lum_ms.append(lum_i)
            if mass >= 1.12:
                R_i = 10 ** (0.66) * np.log(mass) + 0.5
            else:
                R_i = mass
            t_eff_i = (lum_i * L_solar / (4 * np.pi * R_i**2 * R_solar * sigma)) ** (1 / 4)
            t_eff.append(t_eff_i)
        elif t_cool[i] > 0:  # White Dwarf
            if mass < 10:
                num_wd += 1
                mass_wd_i = 0.49 * np.exp(0.095 * mass)
                R_wd_i = C / mass_wd_i**(1 / 3)
                lum_wd_i = 10**(-3) / t_cool[i]**(7 / 5)
                lum_wd.append(lum_wd_i)
                t_eff_wd_i = (lum_wd_i * L_solar / (4 * np.pi * R_wd_i**2 * R_solar * sigma)) ** (1 / 4)
                t_eff_wd.append(t_eff_wd_i)

    # Store results for this cluster age
    all_lum_ms.append(lum_ms)
    all_t_eff.append(t_eff)
    all_lum_wd.append(lum_wd)
    all_t_eff_wd.append(t_eff_wd)

# Plot Main Sequence and White Dwarfs for all ages (normal scale)
plt.figure(figsize=(10, 6))
colors = ['#00008B', '#FF4500', '#008000']  # Colors for ages
for i, age in enumerate(cluster_ages):
    plt.scatter(all_t_eff[i], all_lum_ms[i], s=2, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(all_t_eff_wd[i], all_lum_wd[i], s=2, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Normal Scale)')
plt.gca().invert_xaxis()
plt.legend()
plt.show()

# Plot Main Sequence and White Dwarfs for all ages (log-log scale)
plt.figure(figsize=(10, 6))
for i, age in enumerate(cluster_ages):
    plt.scatter(all_t_eff[i], all_lum_ms[i], s=2, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(all_t_eff_wd[i], all_lum_wd[i], s=2, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Log-Log Scale)')
plt.gca().invert_xaxis()
plt.legend()
plt.show()
