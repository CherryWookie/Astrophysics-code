import numpy as np
import matplotlib.pyplot as plt

np.random.seed(9)

# Create variables
N = 1000 # Number of objects
t_cluster = [8,10,12]
t_born = []
t_ms = []
t_cool = []
dt = 1
alpha = 2.35
min_mass = 0.1
max_mass = 100
mass_ms = []

# Monte Carlo Simulation
for i in range(N):
    t_i = dt * np.random.uniform(0,1)
    t_born.append(t_i)

    # Generate random number for comparison
    rand_y = np.random.uniform(0.1, 100)

    # Generate IMF with random mass vector
    rand_mass = np.random.uniform(0.1, 100)
    mass_i = rand_mass**-alpha

    # Accept Reject
    if mass_i < rand_y:
        mass_ms.append(mass_i)
    else: 
        mass_ms.append(rand_y)

    # Calculate t_ms
    t_msi = 10 * mass_ms[i]**-3.5
    t_ms.append(t_msi)

# Finding t
for years in t_cluster:
    for i in range(N):
        t_cool_i = years - t_born[i] - t_ms[i]
        t_cool.append(t_cool_i)

        if t_cool_i < 0:
            x = 0 # Arbitrary line
            # NO -- MS -- Lum, T_eff
        elif t_cool_i > 0:
            x = 0
            # Yes -- M_ms < 10 -- yes --> White Dwarf --> Lum, T_eff WD. if > 10, no, BH NS




# Plot results
plt.hist(mass_ms, bins=100, log=True, edgecolor='black')
plt.xlabel('Mass (Msun)')
plt.ylabel('Number of Stars')
plt.title('Mass Distribution')
plt.show()

plt.hist(t_born, bins=100, log=True, edgecolor='black')
plt.xlabel('Mass (Msun)')
plt.ylabel('Number of Stars')
plt.title('TBORN')
plt.show()



