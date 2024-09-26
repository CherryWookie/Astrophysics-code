import numpy as np


# Random Number Generator
# This one is bad though because it presents some correlation between numbers, although having a large period

Z = 7 #Seed
a = 2
c = 4
m = 2**31 #full period
iters = 10
output = []

for n in range(iters):
    Z_1 = (a * Z + c) % m
    randu = (65539 * Z) % m #RANDU?
    Z = Z_1
    output.append(randu) # Change append arguement 

print(output)

