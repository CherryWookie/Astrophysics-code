import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


np.random.seed(9)

# Create variables
N = 1000           # Number of objects
# t_cluster = 10      # Cluster age
t_born = []         # Born age
t_ms = []           # Main Sequence Temperature
t_cool = []         # T_cool temperatures
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



# Define cluster ages
cluster_ages = [8, 10, 12]

# Store results for different cluster ages
all_lum_ms = []      # Main sequence luminosities
all_t_eff = []       # Main sequence effective temperatures
all_lum_wd = []      # White dwarf luminosities
all_t_eff_wd = []    # White dwarf effective temperatures

# Simulate for each cluster age
for t_cluster in cluster_ages:

    np.random.seed(9)  # Ensure reproducibility for each age

    # Initialize lists for combined graphs
    mass_ms, t_ms, t_born, t_cool = [], [], [], []
    lum_ms, lum_wd, t_eff, t_eff_wd = [], [], [], []
    num_wd, num_bh = 0, 0

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
                lum_i = 32000 * mass
            elif 2 < mass <= 55:
                lum_i = 1.4 * mass ** 3.5
            elif 0.43 < mass <= 2:
                lum_i = mass ** 4
            elif mass <= 0.43:
                lum_i = 0.23 * mass ** 2.3
            lum_ms.append(lum_i)

            # Main Sequence Radius
            if mass >= 1.12:
                R_i = 10 ** (0.66) * np.log(mass) + 0.5
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

    # Store results for this cluster age
    all_lum_ms.append(lum_ms)
    all_t_eff.append(t_eff)
    all_lum_wd.append(lum_wd)
    all_t_eff_wd.append(t_eff_wd)

    print(f'Number of Black Holes/Neutron Stars for Age: {t_cluster} = {num_bh}')
    print(f'Number of Black White Dwarfs for Age: {t_cluster} = {num_wd}')



    #___________________________________________________________________
    # Individual Plots:
    # Un-comment for individual plots
    #___________________________________________________________________

    # # Create Gaussian Noise Main Sequence
    # lum_noisy = lum_ms * (1 + np.random.normal(0, 0.1, size=len(lum_ms)))
    # temp_noisy = t_eff * (1 + np.random.normal(0, 0.1, size=len(t_eff)))

    # # Create Gaussian Noise White Dwarf
    # lum_wd_noisy = lum_wd * (1 + np.random.normal(0, 0.1, size=len(lum_wd)))
    # temp_wd_noisy = t_eff_wd * (1 + np.random.normal(0, 0.1, size=len(t_eff_wd)))

    # # Debug, Masses from largest to smallest
    # # sorted_arr = sorted(mass_ms, reverse=True)
    # # print("Sorted array from largest to smallest:", sorted_arr)

    # # Print Number of Black Holes and Neutron Stars, as well as White Dwarfs
    # print(f"Number of Black Holes or Neutron Stars: {num_bh}")
    # print(f"Number of White Dwarfs: {num_wd}")

    # # Plot Original Mass 
    # plt.hist(mass_ms, bins=100, log=True, edgecolor='black')
    # plt.xlabel('Mass (Msun)')
    # plt.ylabel('Number of Stars')
    # plt.title('Mass Distribution')
    # plt.show()

    # # Plot t_born
    # plt.hist(t_born, bins=100, log=True, edgecolor='black')
    # plt.xlabel('T_born (Gyrs)')
    # plt.ylabel('Number of Stars')
    # plt.title('TBORN')
    # plt.show()

    # # Plot Main Sequence Luminosity vs T_eff
    # plt.scatter(t_eff, lum_ms, color='#00008B')

    # # Plot White Dwarf Luminosity vs T_eff_wd
    # plt.scatter(t_eff_wd, lum_wd, color='#FFA500')
    # plt.yscale('log')
    # plt.xlabel('T_eff')
    # plt.ylabel('Luminosity')
    # plt.title('Plot of Temperature vs Luminosity')
    # plt.gca().invert_xaxis()
    # plt.show()

    # # Plot White Dwarf / Main Sequence Luminosity vs T_eff_wd LOG-LOG
    # plt.scatter(t_eff, lum_ms, color='#00008B')
    # plt.scatter(t_eff_wd, lum_wd, color='#FFA500')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.xlabel('T_eff')
    # plt.ylabel('Luminosity')
    # plt.title('Plot of Temperature vs Luminosity')
    # plt.gca().invert_xaxis()
    # plt.show()

    # # Plot with 10% Gaussian Noise
    # plt.scatter(temp_noisy, lum_noisy, s=2, color='#00008B')
    # plt.scatter(temp_wd_noisy, lum_wd_noisy, s=2, color='#FFA500')
    # plt.yscale('log')
    # plt.xlabel('T_eff')
    # plt.ylabel('Luminosity')
    # plt.title('Plot of Temperature vs Luminosity')
    # plt.gca().invert_xaxis()
    # plt.show()

    # # Plot with 10% Gaussian Noise LOG-LOG
    # plt.scatter(temp_noisy, lum_noisy, s= 2, color='#00008B')
    # plt.scatter(temp_wd_noisy, lum_wd_noisy, s = 2, color='#FFA500')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.xlabel('T_eff')
    # plt.ylabel('Luminosity')
    # plt.title('Plot of Temperature vs Luminosity')
    # plt.gca().invert_xaxis()
    # plt.show()

#----------------------------------------------------------
# Combined Plots
#----------------------------------------------------------

# Plot Main Sequence and White Dwarfs for all ages (normal scale)
plt.figure(figsize=(10, 6))
colors = ['#00008B', '#FF4500', '#008000']  # Colors for ages
for i, age in enumerate(cluster_ages):
    plt.scatter(all_t_eff[i], all_lum_ms[i], s=20, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(all_t_eff_wd[i], all_lum_wd[i], s=20, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')

plt.yscale('log')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Normal Scale)')
plt.gca().invert_xaxis()
plt.legend()
plt.show()

# Plot Main Sequence and White Dwarfs for all ages (log-log scale)
plt.figure(figsize=(10, 6))
for i, age in enumerate(cluster_ages):
    plt.scatter(all_t_eff[i], all_lum_ms[i], s=10, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(all_t_eff_wd[i], all_lum_wd[i], s=10, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Log-Log Scale)')
plt.gca().invert_xaxis()
plt.legend()
plt.show()

# Plot Gaussian Noise for all ages (normal scale)
plt.figure(figsize=(10, 6))
for i, age in enumerate(cluster_ages):
    # Create Gaussian Noise Main Sequence
    lum_noisy = all_lum_ms[i] * (1 + np.random.normal(0, 0.1, size=len(all_lum_ms[i])))
    temp_noisy = all_t_eff[i] * (1 + np.random.normal(0, 0.1, size=len(all_t_eff[i])))

    # Create Gaussian Noise White Dwarf
    lum_wd_noisy = all_lum_wd[i] * (1 + np.random.normal(0, 0.1, size=len(all_lum_wd[i])))
    temp_wd_noisy = all_t_eff_wd[i] * (1 + np.random.normal(0, 0.1, size=len(all_t_eff_wd[i])))

    # Plot Main Sequence and White Dwarfs with noisy data
    plt.scatter(temp_noisy, lum_noisy, s=10, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(temp_wd_noisy, lum_wd_noisy, s=10, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')

plt.yscale('log')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Normal Scale)\nGaussian Noise 10%')
plt.gca().invert_xaxis()
plt.legend()
plt.show()

# Plot Gaussian Noise for all ages (log-log scale)
plt.figure(figsize=(10, 6))
for i, age in enumerate(cluster_ages):
    # Create Gaussian Noise for Main Sequence
    lum_noisy = all_lum_ms[i] * (1 + np.random.normal(0, 0.1, size=len(all_lum_ms[i])))
    temp_noisy = all_t_eff[i] * (1 + np.random.normal(0, 0.1, size=len(all_t_eff[i])))

    # Create Gaussian Noise for White Dwarf
    lum_wd_noisy = all_lum_wd[i] * (1 + np.random.normal(0, 0.1, size=len(all_lum_wd[i])))
    temp_wd_noisy = all_t_eff_wd[i] * (1 + np.random.normal(0, 0.1, size=len(all_t_eff_wd[i])))

    # Plot Main Sequence and White Dwarfs with gaussian noise
    plt.scatter(temp_noisy, lum_noisy, s=10, color=colors[i], label=f'MS, Age={age} Gyr')
    plt.scatter(temp_wd_noisy, lum_wd_noisy, s=10, color=colors[i], alpha=0.5, label=f'WD, Age={age} Gyr')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('T_eff (K)')
plt.ylabel('Luminosity (L_sun)')
plt.title('Temperature vs Luminosity (Log-Log Scale)\nGaussian Noise 10%')
plt.gca().invert_xaxis()
plt.legend()
plt.show()


#---------------------------------------------
# Adding 2d HR Diagram
#---------------------------------------------


# Combine Main Sequence and White Dwarf Data
t_eff_combined = np.concatenate(all_t_eff)  
lum_combined = np.concatenate(all_lum_ms)   
t_eff_combined_wd = np.concatenate(all_t_eff_wd)  
lum_combined_wd = np.concatenate(all_lum_wd)     

plt.figure(figsize=(10, 6))

# Create a KDE plot (Main Sequence + White Dwarfs)
sns.kdeplot(x=t_eff_combined, y=lum_combined, cmap='inferno', fill=True, thresh=0, levels=30, label='Main Sequence')
sns.kdeplot(x=t_eff_combined_wd, y=lum_combined_wd, cmap='coolwarm', fill=True, thresh=0, levels=30, label='White Dwarfs')

plt.xlabel('Effective Temperature (T_eff) [K]')
plt.ylabel('Luminosity (L/L_sun)')
plt.title('2D KDE Density Map of the HR Diagram')
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_xaxis()
plt.legend()
plt.show()
