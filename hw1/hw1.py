# This is my first assignment for Astrophysics
# Also when getting git to work, I had to change the .json settings to usr/local/bin/git instead of usr/local/bin/git like everything else is
# WTF why was it switched like that. Either way it's working now I suppose
import numpy as np
import math


r,a,b = 2,1,-1


fx1 = b + np.sqrt(r/a)
fx2 = b - np.sqrt(r/a)

N_fp = [0,fx2,fx1]
fx_array = []
words = []

for N in N_fp:
    fx_prime = 3*(N**2) - 4*N*b + b**2 - (r/a)

    print(f'r,a,b: ', r,a,b)
    print(f'N:',N)
    print(f'd/dN: ', fx_prime)

    fx_array.append(fx_prime)
    if fx_prime < 0:
        print("UNSTABLE")
    elif fx_prime > 0:
        print("STABLE")
    else: 
        print("ZERO BOI")
    


