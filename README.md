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
    