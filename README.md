# MethodOfParticularSolutions.jl

This is a package for computing eigenvalues and eigenfunctions of the
Laplacian on planar or spherical domains using the method of
particular solutions in Julia.

It's the implementation of the method described in the article
[Computation of Tight Enclosures for Laplacian
Eigenvalues](https://arxiv.org/abs/2003.08095) and version 0.1.0 were
used to produce the results in it. The article describes the method
and the mathematical background but does not give any code examples.

## Installation
The package is not in the general Julia repositories and does in
addition depend on
[ArbTools.jl](https://github.com/Joel-Dahne/ArbTools.jl) which is not
in the repositories either. You can install both of them through the
package manager.
``` julia
pkg> add https://github.com/Joel-Dahne/ArbTools.jl
pkg> add https://github.com/Joel-Dahne/MethodOfParticularSolutions.jl
```

To see if it works correctly you can run the tests with
``` julia
pkg> test MethodOfParticularSolutions
```
Which should give an output similar to
```
1: Computing eigenvalue for the Spherical triangle with angles (3π/4, 1π/3, 1π/2)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          10    0.25450        [0.000111 +/- 3.27e-7]    [12.40 +/- 7.39e-3]
  16      80           60          18    0.19186          [2.1e-7 +/- 3.12e-9]    [12.4001 +/- 6.67e-5]
radius = 1.827627e-05 < 2.000000e-05

2: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/3, 1π/2)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          17    0.32456          [8.3e-7 +/- 4.43e-9]    [13.7444 +/- 8.93e-5]
  16      80           60          31    0.22173       [4.51e-11 +/- 2.37e-14]    [13.74435521 +/- 6.72e-9]
radius = 3.505609e-09 < 4.000000e-09

3: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/4, 1π/2)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           7    0.21825         [0.00065 +/- 3.73e-6]    [2.1e+1 +/- 0.499]
  16      80           60          15    0.11814      [2.4170e-6 +/- 4.76e-11]    [20.572 +/- 5.09e-4]
radius = 4.815906e-04 < 5.000000e-04

4: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/3, 1π/3)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      70           50          39    0.36021       [4.00e-13 +/- 5.69e-16]    [21.3094076302 +/- 3.38e-11]
  16     120          100          72    0.19384        [2.5e-23 +/- 4.23e-25]    [21.30940763019044525895 +/- 6.29e-21]
radius = 2.804714e-21 < 3.000000e-21

5: Computing eigenvalue for the Spherical triangle with angles (3π/4, 1π/4, 1π/3)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          13    0.33307        [1.689e-5 +/- 2.84e-9]    [24.46 +/- 4.36e-3]
  16      80           60          26    0.20911         [1.3e-9 +/- 2.86e-11]    [24.456914 +/- 3.59e-7]
radius = 1.551502e-07 < 2.000000e-07

6: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/4, 1π/4)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          27    0.18512        [1.17e-9 +/- 6.57e-12]    [49.109945 +/- 4.78e-7]
  16      80           60          49    0.13385       [1.14e-16 +/- 6.81e-19]    [49.1099452632846 +/- 4.03e-14]
radius = 3.026553e-14 < 4.000000e-14

7: Computing eigenvalue for the Spherical triangle with angles (2π/3, 3π/4, 3π/4)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           0    0.17597        [0.038244 +/- 2.98e-7]    [+/- 7.30]
  16      80           60           5    0.13719       [0.0019790 +/- 7.26e-9]    [4e+0 +/- 0.384]
radius = 1.177990e-01 < 1.200000e-01

8: Computing eigenvalue for the Spherical triangle with angles (2π/3, 2π/3, 2π/3)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           2    0.20034        [0.011522 +/- 1.78e-7]    [5e+0 +/- 0.748]
  16      80           60           4    0.15461       [0.0025643 +/- 4.27e-8]    [5e+0 +/- 0.315]
radius = 1.517791e-01 < 2.000000e-01

9: Computing eigenvalue for the Spherical triangle with angles (1π/2, 2π/3, 3π/4)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30    -9223372036854775807    0.11119         [0.10596 +/- 3.00e-6]    nan
  16      80           60           4    0.07624       [0.0010125 +/- 2.20e-8]    [6e+0 +/- 0.386]
radius = 1.407946e-01 < 2.000000e-01

10: Computing eigenvalue for the Spherical triangle with angles (1π/2, 2π/3, 2π/3)
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           1    0.14324        [0.015093 +/- 1.29e-7]    [+/- 8.18]
  16      80           60          12    0.10970      [7.6107e-6 +/- 1.80e-11]    [6.777 +/- 8.70e-4]
radius = 7.607602e-04 < 8.000000e-04

Test Summary:       | Pass  Total
spherical triangles |   20     20
Computing eigenvalue for the L-shaped domain
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           32           3    0.36979          [0.0193 +/- 3.42e-5]    [+/- 10.6]
  16      84           64           8    0.26842     [0.00036790 +/- 7.70e-10]    [9.6 +/- 0.0628]
radius = 2.288494e-02 < 3.000000e-02
Test Summary: | Pass  Total
L-shape       |    2      2
   Testing MethodOfParticularSolutions tests passed
```

## Example
We here show an example of using the package for computing the
fundamental eigenvalue and eigenfunction for the classical example of
the L-shaped domain. The presentation follows the same structure as in
the article mentioned in the beginning which in turn is very similar
to the presentation by Betcke and Trefethen in "Reviving the method of
particular solutions". However the focus here is on the code rather
than the algorithms.

As a first step we load the required packages and set the precision to
use in the computations, in this case 53 bits.

``` julia
using MethodOfParticularSolutions, Nemo, ArbTools, Plots, LaTeXStrings
RR = ArbField(53)
setprecision(BigFloat, 53)
```

We then define the domain that we want to compute the eigenvalues of,
in this case we use the predefined domain for the L-shaped domain.

``` julia
domain = LShape(RR)
```

Next step is to define the particular solution to use. We use a
predefined one which is the same as that used by Betcke and Trefethen.

``` julia
u = LShapeEigenfunction(domain)
```

We can know produce a plot of `sigma(lambda)` similar to Figure 5.2 in
Betcke and Trefethen.

``` julia
N = 15
λs = range(1, stop = 20, length = 200)
σs = [sigma(λ, domain, u, N) for λ in λs]

plot(λs, σs,
     xlims = (0, 20),
     ylims = (0, 1),
     xlabel = L"\lambda",
     ylabel = L"\sigma(\lambda)",
     legend = :none)
```
![Plot of sigma(lambda)](figures/lshape-sigma.png)

We can compute an approximation of the first eigenvalue and
eigenfunction. First we need an interval containing the eigenvalue we
are looking for, since the eigenvalue is approximately given by
9.6397238440219 the interval [9, 10] will do

``` julia
N = 15
interval = setinterval(RR(9), RR(10))
λ = mps(domain, u, interval, N)
```

This produces the, rather poor, enclosure `[9.6 +/- 0.0679]` and also
sets the coefficients of `u` to that of the approximate eigenfunction.
Using a larger value of `N` we can get a better approximation.

``` julia
N = 32
interval = setinterval(RR(9), RR(10))
λ = mps(domain, u, interval, N)
```

This gives us the enclosure `[9.6397 +/- 3.99e-5]` which is slightly
better.

Often times you want to iteratively use higher and higher values of
`N` to get better and better approximations. To achieve this we have
the method `iteratemps`. We can use this to create a figure similar to
Figure 5.3 in Betcke and Trefethen to better see the convergence, we
can plot both the approximate error (computed in the same way as they
do) and the rigorous error given by the radius of the enclosing ball.

``` julia
Ns = 6:4:60
λs = iteratemps(domain, u, interval, Ns,
                optim_prec_final = prec(RR),
                extra_prec = 0,
                show_trace = true)

p = plot(Ns[1:end-1], Float64.(abs.(λs[1:end-1] .- λs[end])),
         xaxis = ("N", 0:10:Ns[end]+9),
         yaxis = ("Error", :log10),
         marker = :cross,
         label = "Approximate error")
plot!(p, Ns, Float64.(radius.(λs)),
      marker = :circle,
      label = "Rigorous error")
```
![Plot of convergence](figures/lshape-convergence.png)

### Higher precision
All of the above computations are done using 53 bits of precision, but
since Arb supports arbitrary precision arithmetic we can go further.
We can use `iteratemps` to compute better and better approximations by
increasing `N` step by step. There are a number of parameters that can
be tuned for this, some of the most important ones are
- values of `N` used, in practice start, step and stop value;
- tolerance used when computing the minimum for each value of `N`;
- precision used in the computation for each value of `N`.

This is better exemplified using a domain for which the convergence is
faster than for the L-shaped domain. We take the spherical triangle
with angles `(2π/3, 1π/3, 1π/3)`, this corresponds to the fourth
triangle in the arXiv paper and we can get the domain, eigenfunction
and an interval containing the first eigenvalue with

``` julia
domain, u, interval = MethodOfParticularSolutions.triangle(4)
```

For this domain taking `N` in steps of 16 works well. We set the
tolerance of the minimization by giving a precision to which to
compute the minimum. In this case we set the precision to be used for
the final value of `N` to 100 and that it should scale linearly with
`N`. Finally we set `extra_prec` to 20 which makes it use a precision
in the computations given by the precision for computing the minimum
plus 20.

``` julia
Ns = 16:16:48
iteratemps(domain, u, interval, Ns,
           optim_prec_final = 250,
           optim_prec_linear = true,
           extra_prec = 20,
           show_trace = true)
```

This takes some time to compute but should give output similar to

```
   N    Prec     Opt prec    Enc prec       Norm                       Maximum    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
  16     104           84          72    0.19384        [2.5e-23 +/- 4.33e-25]    [21.30940763019044525895 +/- 6.29e-21]
  32     187          167         135    0.13469       [1.60e-42 +/- 2.50e-45]    [21.309407630190445258953481441230517778337 +/- 4.17e-40]
  48     270          250         197    0.11049       [2.66e-61 +/- 4.61e-64]    [21.309407630190445258953481441230517778336842577146716613113 +/- 1.96e-58]
3-element Array{arb,1}:
 [21.30940763019044525895 +/- 6.29e-21]
 [21.309407630190445258953481441230517778337 +/- 4.17e-40]
 [21.309407630190445258953481441230517778336842577146716613113 +/- 1.96e-58]
```

Some more information about the options for `iteratemps` are given in
the documentation of the function. However this is not necessarily
complete, the most accurate description you get by reading the actual
code.

## References

Dahne, J., & Salvy, B., Computation of tight enclosures for laplacian
eigenvalues (2020).

Fox, L., P. Henrici, and C. Moler. "Approximations and bounds for
eigenvalues of elliptic operators." SIAM Journal on Numerical Analysis
4.1 (1967): 89-102.

Betcke, Timo, and Lloyd N. Trefethen. "Reviving the method of
particular solutions." SIAM review 47.3 (2005): 469-491.
