# Astrophysics-code
This is the code for my computational astrophysics class, yay!

## Notes

### Random Numbers

Usually random number generators provide a uniform distribution within the interval [0,1]

**Inverse Technique**
- Fundamental law of probability transformation

    - $|p(y)dy| = |p(x)dx|$
- Choosing $p(x)$ according to a uniform distribution and $p(y) = f(y)$
    - $dx/dy = f(y) --> x = F(y)$
    - $x = \int{f(y)dy} = F(y)$
    - $y = F^{-1}(x)$
    - $ x = -e^{-\lambda y}$
        - $ y = -1/\lambda ln|x|$

- Applied to a Gaussian
    - $p(y)dy = e^{-y^{2}/2}/\sqrt{2\pi}  ... $
    - Change variables to y1 and y2

#### Accept-Rejection Method (Monte Carlo)

- Most widely used
- It does not require integration or function inversion

1. Generate random number $x$ between $[a,b]$
2. Generate maximum value $y$ that is larger than any value in the function between $[c,d]$
3. Calculate $p(x):$
    1. if $y > p(x)$ reject it
    2. if $y < p(x)$ accept it

Can be optimized:
- $c = max(f(x))$ otherwise the algorithm slows down
- Choose a known $g(x)$ such that $g(x)>f(x)$, $x$ is selected according to $g(x)$

In order to test for randomness, we can use one of the following:
- *Uniformity test*
    - Values between 0 and 1 should follow uniform distribution
- *Correlation test*
    - We have a set of number $U_1, U_2,...U_n$
    - We make *k*-tuples: ($U_{i+1},U_{i+2},...U_{i+k}$)
    - Evaluate the uniformity distribution on the $(a_1,b_1)(a_2,b_2), ... (a_n,b_n)$
        - Should be able to visualize any uniformity on a plot of the tuples. Noise is better.
    - IMB RANDU Generator:
    $$x_{i+1} = (65539x_i) mod 2^{31}$$
    - RANECU Generator
    - Your own generator (add own constant, i.e. $x_{i+1} = ax_i + b$ mod $c$)


### Uniform Distribution on a Circle

1. We can take into account the number of objects per infinitesimal unit area: $ dn = rd\theta dr$
2. We can also enclose the circle in a box and calculate the distance relative to the center and the radius.
3. If the distribution is uniform the umber of objects per unit volume is constant:
$$ dn = r^2sin\theta dr d\theta d\phi $$

Alternative method is called *Rejection Method*. 

### Astrophysics

Error propagation: The parallax Case

- Parallax
- Units:
    - Astronomical Units (AU)
    - Light year
    - Parsec (pc) = 3.26 light-years
    - Arc sec
        - Distance travelled in one second across earth's surface as an angle(1/3600 of degree)
        - Unit of angular measurement, 3600 arcsec in 1 degree
- Globular Clusters
    - Born at same time
    - Metallicity Z = cte
        - H, He, Z
        - Beginning of "Big Bang" creation, Z was very small. Now the Metallicity has increased over time. Now it is 1%
- Hertzsprung-Russell diagram
    - Blue is hotter, Red is cooler obvi
    - Diagonal is called main sequence track/diagonal
    - Lifetime Main-Sequence

SIMULATION

Create sequence of main-sequence stars and calculate how many will turn into white dwarf stars or something


### Monte Carlo Simulation vs Observed Data

- Observed data is more varied due to unidentifiable objects that could be something other than white dwarfs.
- Luminosity Function


### Simulate Globular Cluster

A globular cluster is a type of stellar aggrupation, formed by ~$10^5 - 10^6$ stars with similar characteristics (formation, metallicity) and gravitationally bound.
 - Globular clusters present great advantages in population synthesis studies:
    - Homogeneous metallicity
    - Short formation period (burst)
    - Bound spacial distribution

#### Simulation Details

Inputs
1. Number of objects $N = 10^4$
2. Cluster age, $t_{cluster} = 8,10,12$ Gyr.
3. Generate for each star it's born time, $t_{born}$, following a burst of star formation of $\delta t = 1.0$ Gyr
$$ t_{born,i} = \delta t * U(0,1)$$

*This should be defined as a vector*

4. Generate mass of the main sequence star, *M_MS* following a IMF. We adopt the IMF from Salpeter (1995) with \alpha = 2.35

$$ \phi(M_{MS}) = (M/M_\omega)^{-\alpha} , with   0.1 < M/M_\omega < 100$$

5. Generate luminosity of the main-sequence stars according to Bressan et al. (1993)