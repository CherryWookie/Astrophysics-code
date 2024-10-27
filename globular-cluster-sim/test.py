import numpy as np
import matplotlib.pyplot as plt

# Set seed for reproducibility
np.random.seed(9)

# Constants and parameters
N = 10**4                           # Number of objects
t_clusters = [8, 10, 12]             # Cluster ages in Gyr
dt = 1.0                             # Star formation burst duration
alpha_salpeter = -2.35               # Slope for Salpeter IMF
sigma = 5.67e-8                      # Stefan-Boltzmann constant (W/m^2 K^4)

# Mass ranges and function parameters
min_mass = 0.1                       # Minimum mass in solar masses
max_mass = 100                       # Maximum mass in solar masses

# Function for Salpeter IMF
def generate_salpeter_mass(N):
    masses = []
    while len(masses) < N:
        rand_mass = np.random.uniform(min_mass, max_mass)
        phi_mass = rand_mass ** alpha_salpeter
        rand_y = np.random.uniform(0, max_mass**(alpha_salpeter + 1))
        if rand_y < phi_mass:
            masses.append(rand_mass)
    return masses

# Main sequence luminosity based on Salaris & Cassisi (2005)
def calculate_luminosity(mass):
    if mass <= 0.43:
        return 0.23 * mass ** 2.3
    elif 0.43 < mass <= 2:
        return mass ** 4
    elif 2 < mass <= 55:
        return 1.4 * mass ** 3.5
    else:
        return 32000 * mass

# Main sequence lifetime based on Iben & Laughlin (1989)
def calculate_t_ms(mass):
    return 10 * mass ** -3.5  # Gyr

# White dwarf mass from initial-final mass relation
def calculate_wd_mass(mass):
    return 0.49 * np.exp(0.095 * mass)

# Calculate cooling age
def calculate_cooling_age(t_cluster, t_born, t_ms):
    return t_cluster - t_born - t_ms

# Calculate luminosity using Mestel's cooling model
def calculate_luminosity_cooling(t_cool):
    return 10 ** (-5 / 7 * np.log10(t_cool) + 3)

# Effective temperature from Stefan-Boltzmann law
def calculate_t_eff(lum, radius):
    return (lum / (4 * np.pi * radius**2 * sigma)) ** (1/4)

# Main-sequence radius
def calculate_radius(mass):
    if mass >= 1.12:
        return 10 ** (0.66 * np.log10(mass)) + 0.05
    else:
        return mass

# Monte Carlo simulation function
def run_simulation(t_cluster, mass_function='salpeter'):
    # Initialize arrays
    masses, lum, t_born, t_ms, t_cool, wd_lum, t_eff, R, wd_masses = [], [], [], [], [], [], [], [], []
    num_wd, num_ns_bh = 0, 0  # Counters for WDs and neutron stars/black holes

    # Generate initial masses based on chosen IMF
    mass_ms = generate_salpeter_mass(N)

    for mass in mass_ms:
        # Random age at birth within [0, dt]
        t_i = dt * np.random.uniform(0, 1)
        t_born.append(t_i)

        # Main sequence lifetime
        t_msi = calculate_t_ms(mass)
        t_ms.append(t_msi)

        # Calculate MS luminosity and radius
        lum_i = calculate_luminosity(mass)
        lum.append(lum_i)

        R_i = calculate_radius(mass)
        R.append(R_i)

        # Calculate cooling time
        t_cool_i = calculate_cooling_age(t_cluster, t_i, t_msi)
        t_cool.append(t_cool_i)

        # Post-MS classification and cooling
        if t_cool_i > 0:
            if mass < 10:
                # White Dwarf
                wd_mass = calculate_wd_mass(mass)
                wd_masses.append(wd_mass)
                wd_lumi = calculate_luminosity_cooling(t_cool_i)
                wd_lum.append(wd_lumi)
                wd_temp = calculate_t_eff(wd_lumi, R_i)
                t_eff.append(wd_temp)
                num_wd += 1
            else:
                # Neutron star / black hole
                num_ns_bh += 1
        else:
            # Main sequence star properties
            t_eff_i = calculate_t_eff(lum_i, R_i)
            t_eff.append(t_eff_i)

    return {
        "lum": lum, "t_eff": t_eff, "R": R, "wd_lum": wd_lum, "wd_masses": wd_masses,
        "num_wd": num_wd, "num_ns_bh": num_ns_bh
    }

# Run simulation for each cluster age
for t_cluster in t_clusters:
    print(f"Running simulation for t_cluster = {t_cluster} Gyr")
    results = run_simulation(t_cluster, mass_function='salpeter')
    print(f"Total WDs: {results['num_wd']}, NS/BH: {results['num_ns_bh']}")

    # Plotting: Check that we have data to plot
    if results['t_eff'] and results['lum']:
        plt.figure()
        plt.plot(results['t_eff'], results['lum'], color='blue', label='Main Sequence')
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Effective Temperature (K)')
        plt.ylabel('Luminosity (L/Lsun)')
        plt.title(f'L-T Plot for Cluster Age {t_cluster} Gyr')
        plt.legend()
        plt.gca().invert_xaxis()
        plt.show()
    else:
        print(f"No data available for plotting for t_cluster = {t_cluster} Gyr.")
