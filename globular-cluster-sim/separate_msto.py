import numpy as np
import matplotlib.pyplot as plt

# This code is to find the MSTO 

# Constants
np.random.seed(9)
N = 1000           # Number of stars
dt = 1             # Timestep
alpha = 2.35       # IMF constant
max_mass = 100     # Maximum star mass
sigma = 5.67e-8    # Stefan-Boltzmann constant (W/m^2 K^4)
C = 0.0101         # White Dwarf constant
R_solar = 6.969e8  # Solar radius (m)
L_solar = 3.827e26 # Solar luminosity (W)

# Cluster ages and minimum mass values
cluster_ages = [8, 10, 12]
min_mass_values = [0.1, 0.3, 0.5] # Chosen as solar masses

# Main code but split into definitions to simplify
def generate_masses(min_mass, max_mass, alpha, N):
    masses = []
    while len(masses) < N:
        rand_y = np.random.uniform(0, 1)
        rand_mass = np.random.uniform(min_mass, max_mass)
        mass_i = (rand_mass) ** (-alpha)
        if rand_y < mass_i:
            masses.append(rand_mass)
    return masses

def main_sequence_luminosity(mass):
    if mass > 55:
        return 32000 * mass
    elif 2 < mass <= 55:
        return 1.4 * mass ** 3.5
    elif 0.43 < mass <= 2:
        return mass ** 4
    elif mass <= 0.43:
        return 0.23 * mass ** 2.3

def main_sequence_radius(mass):
    if mass >= 1.12:
        return 10 ** 0.66 * np.log(mass) + 0.5
    else:
        return mass

def effective_temperature(lum, radius):
    return (lum * L_solar / (4 * np.pi * radius ** 2 * R_solar ** 2 * sigma)) ** (1 / 4)

def calculate_msto(cluster_age, masses, t_ms):
    # MSTO is the most massive star that is still on the main sequence
    for i, age in enumerate(t_ms):
        if cluster_age - age <= 0:
            return i, masses[i]
    return len(masses) - 1, masses[-1]  # Default to least massive star if no MSTO

# Main simulation
for min_mass in min_mass_values:
    print(f"Analyzing for min_mass = {min_mass} M_sun")
    
    all_lum_ms, all_t_eff, all_msto = [], [], []
    colors = ['#00008B', '#FF4500', '#008000']

    for age_idx, t_cluster in enumerate(cluster_ages):
        # Generate initial star parameters
        masses = generate_masses(min_mass, max_mass, alpha, N)
        t_ms = [10 * mass ** -3.5 for mass in masses]
        t_born = [dt * np.random.uniform(0, 1) for _ in range(N)]
        t_cool = [t_cluster - t_born[i] - t_ms[i] for i in range(N)]

        # Separate stars into MS and WD
        lum_ms, t_eff_ms = [], []
        for i, mass in enumerate(masses):
            if t_cool[i] <= 0:
                lum = main_sequence_luminosity(mass)
                radius = main_sequence_radius(mass)
                t_eff = effective_temperature(lum, radius)
                lum_ms.append(lum)
                t_eff_ms.append(t_eff)

        # Determine MSTO
        msto_index, msto_mass = calculate_msto(t_cluster, masses, t_ms)
        msto_lum = main_sequence_luminosity(msto_mass)
        msto_radius = main_sequence_radius(msto_mass)
        msto_t_eff = effective_temperature(msto_lum, msto_radius)
        all_msto.append((msto_mass, msto_t_eff, msto_lum))

        # Append data
        all_lum_ms.append(lum_ms)
        all_t_eff.append(t_eff_ms)

        # Print MSTO details
        print(f"Age: {t_cluster} Gyr, MSTO Mass: {msto_mass:.2f} M_sun, "
              f"T_eff: {msto_t_eff:.2f} K, Luminosity: {msto_lum:.2f} L_sun")

    # Plot results
    plt.figure(figsize=(10, 6))
    for i, age in enumerate(cluster_ages):
        plt.scatter(all_t_eff[i], all_lum_ms[i], s=10, color=colors[i], label=f'Age = {age} Gyr')
        msto_mass, msto_t_eff, msto_lum = all_msto[i]
        plt.scatter(msto_t_eff, msto_lum, s=100, edgecolor='black', color=colors[i], marker='x', label=f'MSTO (Age={age} Gyr)')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('T_eff (K)')
    plt.ylabel('Luminosity (L_sun)')
    plt.title(f'Temperature vs Luminosity (Log-Log Scale) for min_mass = {min_mass} M_sun')
    plt.gca().invert_xaxis()
    plt.legend()
    plt.show()
