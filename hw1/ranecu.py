import numpy as np

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
rng = RandomNumberGenerator(seed1=12345, seed2=67890)

# Generate a random number
random_number = rng.ran()
print(f'RANECU NUMBER: ', random_number)
