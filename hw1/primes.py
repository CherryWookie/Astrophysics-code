import numpy as np

# This is a program to compute prime numbers up to a certain limit
# Adjust the 'limit' variable depending on how many primes you want

def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def find_primes_up_to(limit):
    primes = []
    for num in range(2, limit + 1):
        if is_prime(num):
            primes.append(num)
    return primes

# Find and print all prime numbers up to 100
limit = 1000
primes_up_to_100 = find_primes_up_to(limit)
print(primes_up_to_100)

# This is something new
