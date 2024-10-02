class RANGenerator:
    def __init__(self, seed1, seed2):
        self.ISEED1 = seed1
        self.ISEED2 = seed2
        self.USCALE = 1.0 / 2147483563

    def ran(self):
        def mod32(val):
            if val > 2147483647:
                val -= 2147483647
            elif val < 0:
                val += 2147483647
            return val

        I1 = self.ISEED1 // 53668
        self.ISEED1 = 40014 * (self.ISEED1 - I1 * 53668) - I1 * 12211
        if self.ISEED1 < 0:
            self.ISEED1 += 2147483563

        I2 = self.ISEED2 // 52774
        self.ISEED2 = 40692 * (self.ISEED2 - I2 * 52774) - I2 * 3791
        if self.ISEED2 < 0:
            self.ISEED2 += 2147483399

        IZ = self.ISEED1 - self.ISEED2
        if IZ < 1:
            IZ += 2147483562

        return IZ * self.USCALE


# Target numbers to match
target_numbers = [
    0.919494450092316, 0.311697453260422, 0.214007109403610,
    0.218347966670990, 0.281165540218353, 0.892118871212006,
    0.358745545148849, 0.422309160232544, 0.139809057116508,
    0.616842269897461
]

# Brute-force search over a range of seeds
def find_seeds(target_numbers):
    for seed1 in range(100000, 500000):  # Searching seed1 in a reasonable range
        for seed2 in range(100000, 500000):  # Searching seed2 in a reasonable range
            rng = RANGenerator(seed1, seed2)
            first_random_number = rng.ran()
            
            # Check if the first generated random number matches the target
            if abs(first_random_number - target_numbers[0]) < 1e-9:
                # Check if subsequent numbers match
                match = True
                for i in range(1, len(target_numbers)):
                    random_number = rng.ran()
                    if abs(random_number - target_numbers[i]) >= 1e-9:
                        match = False
                        break
                
                if match:
                    print(f"Found matching seeds! seed1={seed1}, seed2={seed2}")
                    return seed1, seed2

    print("No matching seeds found in the searched range.")
    return None, None

# Search for the seeds
seed1, seed2 = find_seeds(target_numbers)
