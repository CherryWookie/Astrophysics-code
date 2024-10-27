import numpy as np
import matplotlib.pyplot as plt

# Random Number Generator
# This one is bad though because it presents some correlation between numbers, although having a large period

Z_0 = 434 # Initial Seed
Z_randu = 47934
a = 2 
c = 4
m = 237725 # Full Period for RANDOM NUMBER
n = 2**31 # Full Period for RANDU IBM
iters = 10
output_rand = []
output_randu = []

# Random Number with Z seed and N iterations

for i in range(iters):
    # Random number calculation
    Z_1 = (a * Z_0 + c) % m

    # RANDU IBM number calculation
    randu = (65539 * Z_randu) % n #RANDU?

    Z_0 = Z_1
    Z_randu = randu

    # Append to respective lists
    output_rand.append(Z_1)
    output_randu.append(randu)


print(output_rand)

# Normalize randu output list
output_randu = [x / 2**31 for x in output_randu]

# Normalize rand output list
output_rand = [x / 2**31 for x in output_rand]

# Create pairs (a,b) for 2-d scatter plot
doubles = [(output_randu[i],output_randu[i+1]) for i in range(0, len(output_randu), 3)]

a_vals = [point[0] for point in doubles]
b_vals = [point[1] for point in doubles]

# Create tuples (x,y,z) for randu output list from (a,b,c,d...n) --> [(a,b,c)(d,e,f)...(n,n,n)]
tuples = [(output_randu[i],output_randu[i+1],output_randu[i+2]) for i in range(0, len(output_randu), 4)]

# Create x,y,z values for points by getting index [0,1,2] for each point (x,y,z)
x_vals = [point[0] for point in tuples]
y_vals = [point[1] for point in tuples]
z_vals = [point[2] for point in tuples]

# Create x, y, z for random number generator (MINE :))
rn_tuples = [(output_rand[i], output_rand[i+1], output_rand[i+2]) for i in range(0, len(output_rand) - 2, 3)]

x_rn = [point[0] for point in rn_tuples]
y_rn = [point[1] for point in rn_tuples]
z_rn = [point[2] for point in rn_tuples]

# 3D Plot Random Number Generator (Custom Generator)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_rn, y_rn, z_rn, c='r', label='Random Number Generator')
ax.set_title('RANDom Number Generator')

# Plot points for randum number generator
# plt.scatter(pointsa_rn,pointsb_rn)
# plt.scatter(a_vals,b_vals)

# Create a 3D scatter plot RANDU
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_vals, y_vals, z_vals)
ax.set_title('RANDU Number Generator')

# Adding labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('RANDU Num Generator')

# Show the plot
plt.show()