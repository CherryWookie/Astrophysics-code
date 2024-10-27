import numpy as np
import matplotlib.pyplot as plt

np.random.seed(9)

# Create variables
N = 1000 # Number of objects
t_cluster = 10
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
sigma = 5.67e-8 # W/m^2 K^4

# Monte Carlo Simulation
while len(mass_ms) < N:

    # Generate random number for comparison
    rand_y = np.random.uniform(min_mass, max_mass**(-alpha))

    # Generate IMF with random mass vector between min-max mass
    rand_mass = np.random.uniform(min_mass, max_mass)
    mass_i = rand_mass**(-alpha)
    

    # Accept Reject
    if rand_y < mass_i:
        # ACCEPT
        mass_ms.append(mass_i)
    else: 
        # REJECT
        continue


# Create T_ms and T_i arrays
for i in range(N):
    # Calculate t_ms
    t_msi = 10 * mass_ms[i] ** -3.5
    t_ms.append(t_msi)

    # Calculate t_born
    t_i = dt * np.random.uniform(0, 1)
    t_born.append(t_i)

    # Calculate t_cool
    t_cool_i = t_cluster - t_born[i] - t_ms[i]
    t_cool.append(t_cool_i)

# Finding luminosity and temperature, etc.
for i in range(N):

    mass = mass_ms[i]
    # print(mass)

    # T_cool conditions to calculate Luminosity
    if t_cool[i] < 0:
        # Calculate Luminosity based on mass ranges
        if mass >= 20:
            lum_i = 81 * mass ** 2.14
        elif 2 < mass <= 20:
            lum_i = 1.78 * mass ** 3.5
        elif mass <= 2:
            lum_i = 0.75 * mass ** 4.8
        lum.append(lum_i)

        # Calculate Radius R based on mass
        if mass >= 1.2:
            R_i = 10 ** (0.66 * np.log(mass)) + 0.5
        else:
            R_i = mass
        R.append(R_i)

        # Calculate Effective Temperature T_eff
        t_eff_i = (lum_i / (4 * np.pi * R_i**2 * sigma)) ** (1 / 4)
        t_eff.append(t_eff_i)

        # Debug: Print values in the loop
        # print(f"mass[{i}] = {mass}, t_cool[{i}] = {t_cool[i]}, lum[{i}] = {lum_i}, R[{i}] = {R_i}, T_eff[{i}] = {t_eff_i}")
    else:
        print(f"mass[{i}] = {mass} did not meet t_cool condition, t_cool[{i}] = {t_cool[i]}")

# # Debug: Final output values
# print("Final Luminosity:", lum)
# print("Final Radius:", R)
# print("Final T_eff:", t_eff)



    # elif t_cool_i > 0:
    #     if mass_ms[i] < 10:
    #         num_wd += 1

    #     elif mass_ms[i] > 10:
    #         num_bh += 1

    #     Yes -- M_ms < 10 -- yes --> White Dwarf --> Lum, T_eff WD. if > 10, no, BH (Black hole) NS (Neutron star)

# print(f"Luminosity: {lum}")
# print(f"Radius: {R}")
# print(f"MASSES: {mass_ms}")

sorted_arr = sorted(mass_ms, reverse=True)

# print("Sorted array from largest to smallest:", sorted_arr)
# print(f"Number of Black Holes or Neutron Stars: {num_bh}")
# print(f"Number of White Dwarfs: {num_wd}")
# # print(len(t_ms))

# Plot Mass 
plt.hist(mass_ms, bins=100, log=True, edgecolor='black')
plt.xlabel('Mass (Msun)')
plt.ylabel('Number of Stars')
plt.title('Mass Distribution')
plt.show()

# Plot t_born
plt.hist(t_born, bins=100, log=True, edgecolor='black')
plt.xlabel('T_born (Gyrs)')
plt.ylabel('Number of Stars')
plt.title('TBORN')
plt.show()

# Plot Luminosity vs T_eff
plt.scatter(t_eff, lum)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('T_eff')
plt.ylabel('Luminosity')
plt.title('Plot of Temperature vs Luminosity')
plt.gca().invert_yaxis()
plt.show()


