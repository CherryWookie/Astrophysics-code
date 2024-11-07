import numpy as np
import matplotlib.pyplot as plt

np.random.seed(9)

# Create variables
N = 1000           # Number of objects
t_cluster = 10      # Cluster age
t_born = []         # Born age
t_ms = []           # Main Sequence Temperature
t_cool = []         # T_cool temperature
dt = 1              # Timestep
alpha = 2.35        # IMF constant
min_mass = 0.1      # Minimum Star Mass
max_mass = 100      # Maximum Star Mass
mass_ms = []        # Main Sequence Mass Array 
mass_wd = []        # White Dwarf Mass Array
num_wd = 0          # White Dwarf Count
num_bh = 0          # Black Holes and Neutron Star count
lum_ms = []         # Main Sequence Luminosity Array
lum_wd = []         # White Dwarf Luminosity Array
t_eff = []          # Main Sequence Effective Temperature 
t_eff_wd = []       # White Dwarf Effective Temperature
R_ms = []           # Main Sequence Radius
R_wd = []           # White Dwarf Radius
sigma = 5.67e-8     # W/m^2 K^4
C = 0.0101          # Constant
R_solar = 6.969e8   # Solar Radius  
L_solar = 3.827e26  # Solar Luminosity

# Monte Carlo Simulation
while len(mass_ms) < N:

    # Generate random number for comparison between 1 and 0
    rand_y = np.random.uniform(0,1)

    # Generate IMF with random mass vector between min-max mass
    rand_mass = np.random.uniform(min_mass, max_mass)
    mass_i = (rand_mass)**(-alpha) # Use mass_i as the probability function of acceptance 


    # Accept Reject Method
    if rand_y < mass_i:
        # ACCEPT
        mass_ms.append(rand_mass) # Accept the random mass in between min and max, not mass_i
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

# Finding luminosity and temperature, etc. for MS and WD
for i in range(N):

    mass = mass_ms[i]

    # MAIN SEQUENCE CONDITIONS
    if t_cool[i] <= 0:
        # Calculate Luminosity based on mass ranges
        if mass > 55:
            lum_i = 3200 * mass
        elif 2 < mass <= 55:
            lum_i = 1.4 * mass ** 3.5
        elif 0.43 < mass <= 2:
            lum_i = mass ** 4
        elif mass <= 0.43:
            lum_i = 0.23 * mass ** 2.3
        lum_ms.append(lum_i)

        # Main Sequence Radius
        if mass >= 1.2:
            R_i = 10 ** (0.66 * np.log(mass)) + 0.5
        else:
            R_i = mass
        R_ms.append(R_i)

        # Main Sequence Effective Temperature
        t_eff_i = (lum_i * L_solar/ (4 * np.pi * R_i**2 * R_solar * sigma)) ** (1 / 4)
        t_eff.append(t_eff_i)

        # print(f"mass[{i}] = {mass}, t_cool[{i}] = {t_cool[i]}, lum[{i}] = {lum_i}, R[{i}] = {R_i}, T_eff[{i}] = {t_eff_i}")
 
 
    # WHITE DWARF CONDITIONS
    elif t_cool_i > 0:
        if mass_ms[i] < 10:
            num_wd += 1 # White Dwarf

            # White Dwarf Mass
            mass_wd_i = 0.49 * np.exp(0.095 * mass_ms[i])
            mass_wd.append(mass_wd_i)
            
            # White Dwarf Radius
            R_wd_i = C / mass_wd_i**(1/3)
            R_wd.append(R_wd_i)

            # White Dwarf Luminosity
            lum_wd_i = 10**(-3) / t_cool[i]**(7/5)
            lum_wd.append(lum_wd_i)

            # White Dwarf T_eff
            t_eff_wd_i = (lum_wd_i * L_solar / (4 * np.pi * R_wd_i**2 * R_solar * sigma)) ** (1 / 4)
            t_eff_wd.append(t_eff_wd_i)

        # Black Hole or Neutron Star
        elif mass_ms[i] > 10:
            num_bh += 1 # Black Hole or Neutron Star

#___________________________________________________________________
# PLOTS:
#___________________________________________________________________

# Create Gaussian Noise Main Sequence
lum_noisy = lum_ms * (1 + np.random.normal(0, 0.1, size=len(lum_ms)))
temp_noisy = t_eff * (1 + np.random.normal(0, 0.1, size=len(t_eff)))

# Create Gaussian Noise White Dwarf
lum_wd_noisy = lum_wd * (1 + np.random.normal(0, 0.1, size=len(lum_wd)))
temp_wd_noisy = t_eff_wd * (1 + np.random.normal(0, 0.1, size=len(t_eff_wd)))

# Debug, Masses from largest to smallest
# sorted_arr = sorted(mass_ms, reverse=True)
# print("Sorted array from largest to smallest:", sorted_arr)

# Print Number of Black Holes and Neutron Stars, as well as White Dwarfs
print(f"Number of Black Holes or Neutron Stars: {num_bh}")
print(f"Number of White Dwarfs: {num_wd}")

# Plot Original Mass 
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

# Plot Main Sequence Luminosity vs T_eff
plt.scatter(t_eff, lum_ms, color='#00008B')

# Plot White Dwarf Luminosity vs T_eff_wd
plt.scatter(t_eff_wd, lum_wd, color='#FFA500')
plt.yscale('log')
plt.xlabel('T_eff')
plt.ylabel('Luminosity')
plt.title('Plot of Temperature vs Luminosity')
plt.gca().invert_xaxis()
plt.show()

# Plot White Dwarf / Main Sequence Luminosity vs T_eff_wd LOG-LOG
plt.scatter(t_eff, lum_ms, color='#00008B')
plt.scatter(t_eff_wd, lum_wd, color='#FFA500')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('T_eff')
plt.ylabel('Luminosity')
plt.title('Plot of Temperature vs Luminosity')
plt.gca().invert_xaxis()
plt.show()

# Plot with 10% Gaussian Noise
plt.scatter(temp_noisy, lum_noisy, s=2, color='#00008B')
plt.scatter(temp_wd_noisy, lum_wd_noisy, s=2, color='#FFA500')
plt.yscale('log')
plt.xlabel('T_eff')
plt.ylabel('Luminosity')
plt.title('Plot of Temperature vs Luminosity')
plt.gca().invert_xaxis()
plt.show()

# Plot with 10% Gaussian Noise LOG-LOG
plt.scatter(temp_noisy, lum_noisy, s= 2, color='#00008B')
plt.scatter(temp_wd_noisy, lum_wd_noisy, s = 2, color='#FFA500')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('T_eff')
plt.ylabel('Luminosity')
plt.title('Plot of Temperature vs Luminosity')
plt.gca().invert_xaxis()
plt.show()