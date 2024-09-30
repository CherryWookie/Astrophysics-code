import numpy as np


# Random Number Generator
# This one is bad though because it presents some correlation between numbers, although having a large period

Z_0 = 7 # Initial Seed
a = 2
c = 4
m = 23 # Full Period for RANDOM NUMBER
n = 2**31 # Full Period for RANDU IBM
iters = 10
output_rand = []
output_randu = []

# Random Number with Z seed and N iterations

for i in range(iters):
    Z_1 = (a * Z_0 + c) % m
    randu = (65539 * Z_0) % n #RANDU?
    Z_0 = Z_1
    output_rand.append(Z_1)
    output_randu.append(randu) # Change append arguement 

print(f'RANDOM NUMBER: ', output_rand)
print(f'RANDU IBM: ', output_randu)


