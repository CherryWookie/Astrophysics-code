import numpy as np
import matplotlib.pyplot as plt

# This code is a translation of the Fortran code for a RANECU Random Number generator with 2 seeds
# Michael Sell

iters = 2000 # Number of random numbers outputed

class RandomNumberGenerator:
    def __init__(self, seed1, seed2):
        self.ISEED1 = seed1
        self.ISEED2 = seed2
        self.USCALE = 1.0 / 2147483563.0  # Equivalent of 1.0D0/2.147483563D9

    def ran(self):
        # First seed calculations (similar to I1 and ISEED1 update in Fortran)
        I1 = self.ISEED1 // 53668
        self.ISEED1 = 40014 * (self.ISEED1 - I1 * 53668) - I1 * 12211
        if self.ISEED1 < 0:
            self.ISEED1 += 2147483563

        # Second seed calculations (similar to I2 and ISEED2 update in Fortran)
        I2 = self.ISEED2 // 52774
        self.ISEED2 = 40692 * (self.ISEED2 - I2 * 52774) - I2 * 3791
        if self.ISEED2 < 0:
            self.ISEED2 += 2147483399

        # Calculate the difference and ensure IZ is positive
        IZ = self.ISEED1 - self.ISEED2
        if IZ < 1:
            IZ += 2147483562

        # Generate the random number using the scaling factor
        return IZ * self.USCALE

# Example usage
rng = RandomNumberGenerator(seed1=348713, seed2=738131)

x = []
y = []
ranecu = []
# Generate a random number
for i in range(iters):      
    random_number = rng.ran()
    print(f'RANECU NUMBER: ', random_number)
    ranecu.append(random_number)
    # if i % 2 == 0:
    #     x.append(random_number)
    # else:
    #     y.append(random_number)

# Create tuples for 3-D plot
tuples = [(ranecu[i],ranecu[i+1],ranecu[i+2]) for i in range(0, len(ranecu), 4)]

x_vals = [point[0] for point in tuples]
y_vals = [point[1] for point in tuples]
z_vals = [point[2] for point in tuples]

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_vals, y_vals, z_vals)

# Plot points for randum number generator
# plt.scatter(x,y)

# Adding labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Scatter Plot of Points RANECU')

# Show the plot
plt.show()
