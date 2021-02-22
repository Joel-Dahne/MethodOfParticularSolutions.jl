# MethodOfParticularSolutions.jl

This is a package for computing eigenvalues and eigenfunctions of the
Laplacian on planar or spherical domains using the method of
particular solutions in Julia.

It's the implementation of the method described in the article
[Computation of Tight Enclosures for Laplacian
Eigenvalues](https://arxiv.org/abs/2003.08095) and version 0.1.0 were
used to produce the results in it. The article describes the method
and the mathematical background but does not give any code examples.

The code has since then been reworked for a new project where an
eigenvalue problem is solved in the plane. Some parts of the package
have changed quite a bit whereas other parts are left in their old
state. If you try to work with the code you are likely to find many
inconsistencies.

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
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          10    0.25450        [0.000437 +/- 8.68e-7]    [12.40 +/- 7.38e-3]
  16      80           60          18    0.19186         [1.08e-6 +/- 6.21e-9]    [12.4001 +/- 6.66e-5]
radius = 1.823349e-05 < 2.000000e-05

2: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/3, 1π/2)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          17    0.32456         [2.57e-6 +/- 7.79e-9]    [13.7444 +/- 8.92e-5]
  16      80           60          31    0.22173     [2.0347e-10 +/- 8.26e-15]    [13.74435521 +/- 6.72e-9]
radius = 3.505115e-09 < 4.000000e-09

3: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/4, 1π/2)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           7    0.21825         [0.00297 +/- 8.68e-6]    [2.1e+1 +/- 0.499]
  16      80           60          15    0.11814      [2.0459e-5 +/- 3.39e-10]    [20.572 +/- 5.09e-4]
radius = 4.815291e-04 < 5.000000e-04

4: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/3, 1π/3)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      70           50          37    0.36021       [3.82e-12 +/- 6.11e-15]    [21.3094076302 +/- 9.30e-11]
  16     120          100          68    0.19384       [1.31e-21 +/- 8.67e-24]    [21.3094076301904452590 +/- 7.52e-20]
radius = 2.867225e-20 < 3.000000e-20

5: Computing eigenvalue for the Spherical triangle with angles (3π/4, 1π/4, 1π/3)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          13    0.33307      [5.0713e-5 +/- 5.58e-10]    [24.46 +/- 4.36e-3]
  16      80           60          26    0.20911         [6.1e-9 +/- 9.92e-11]    [24.456914 +/- 3.59e-7]
radius = 1.551478e-07 < 2.000000e-07

6: Computing eigenvalue for the Spherical triangle with angles (2π/3, 1π/4, 1π/4)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30          25    0.18512        [2.43e-8 +/- 5.37e-11]    [49.10995 +/- 5.62e-6]
  16      80           60          44    0.13385        [3.1e-14 +/- 5.14e-16]    [49.10994526328 +/- 5.71e-12]
radius = 1.095922e-12 < 2.000000e-12

7: Computing eigenvalue for the Spherical triangle with angles (2π/3, 3π/4, 3π/4)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           0    0.17598          [0.2173 +/- 7.06e-5]    [+/- 7.30]
  16      80           60           5    0.13720        [0.014425 +/- 4.60e-7]    [4e+0 +/- 0.384]
radius = 1.177888e-01 < 1.200000e-01

8: Computing eigenvalue for the Spherical triangle with angles (2π/3, 2π/3, 2π/3)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           2    0.20034         [0.05751 +/- 8.58e-6]    [5e+0 +/- 0.748]
  16      80           60           4    0.15461        [0.016585 +/- 2.04e-8]    [5e+0 +/- 0.315]
radius = 1.517747e-01 < 2.000000e-01

9: Computing eigenvalue for the Spherical triangle with angles (1π/2, 2π/3, 3π/4)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30    -9223372036854775807    0.11026           [0.915 +/- 2.33e-4]    [+/- inf]
  16      80           60           5    0.07624        [0.011774 +/- 3.74e-7]    [6e+0 +/- 0.370]
radius = 1.247642e-01 < 2.000000e-01

10: Computing eigenvalue for the Spherical triangle with angles (1π/2, 2π/3, 2π/3)
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           30           2    0.14319          [0.0883 +/- 4.50e-5]    [+/- 7.91]
  16      80           60          12    0.10970      [7.3227e-5 +/- 4.11e-10]    [6.78 +/- 3.70e-3]
radius = 8.029671e-04 < 9.000000e-04

Test Summary:       | Pass  Total
spherical triangles |   20     20
Computing eigenvalue for the L-shaped domain
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
   8      53           32           3    0.36979          [0.0522 +/- 5.85e-5]    [+/- 10.6]
  16      84           64           8    0.26842       [0.0013706 +/- 2.18e-8]    [9.6 +/- 0.0628]
radius = 2.288490e-02 < 3.000000e-02
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
λ = iteratemps(domain, u, interval, N:N, optim_prec_final = prec(RR), extra_prec = 0)[1]
```

This produces the, rather poor, enclosure `[9.6 +/- 0.0679]` and also
sets the coefficients of `u` to that of the approximate eigenfunction.
Using a larger value of `N` we can get a better approximation.

``` julia
N = 32
interval = setinterval(RR(9), RR(10))
λ = iteratemps(domain, u, interval, N:N, optim_prec_final = prec(RR), extra_prec = 0)[1]
```

This gives us the enclosure `[9.6397 +/- 3.99e-5]` which is slightly
better.

Often times you want to iteratively use higher and higher values of
`N` to get better and better approximations. To achieve this you can
pass `iteratemps` a range of `N` values to use. We can use this to
create a figure similar to Figure 5.3 in Betcke and Trefethen to
better see the convergence, we can plot both the approximate error
(computed in the same way as they do) and the rigorous error given by
the radius of the enclosing ball.

``` julia
Ns = 6:4:60
λs = iteratemps(domain, u, interval, Ns,
                optim_prec_final = prec(RR),
                extra_prec = 0,
                show_trace = true)[1]

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
           show_trace = true)[1]
```

This takes some time to compute but should give output similar to

```
   N    Prec     Opt prec    Enc prec       Norm                  Maximum/Norm    Enclosure
----    ----     --------    --------    -------    --------------------------    ---------
  16     104           84          68    0.19384       [1.31e-21 +/- 8.75e-24]    [21.3094076301904452590 +/- 7.52e-20]
  32     187          167         129    0.13469       [1.01e-39 +/- 5.89e-42]    [21.3094076301904452589534814412305177783 +/- 5.90e-38]
  48     270          250         188    0.11049       [1.65e-57 +/- 9.46e-60]    [21.3094076301904452589534814412305177783368425771467166131 +/- 4.94e-56]
3-element Array{arb,1}:
 [21.3094076301904452590 +/- 7.52e-20]
 [21.3094076301904452589534814412305177783 +/- 5.90e-38]
 [21.3094076301904452589534814412305177783368425771467166131 +/- 4.94e-56]
```

Some more information about the options for `iteratemps` are given in
the documentation of the function.

## References

Dahne, J., & Salvy, B., Computation of tight enclosures for laplacian
eigenvalues, SIAM} Journal on Scientific Computing, 42(5), 3210–3232
(2020). http://dx.doi.org/10.1137/20m1326520

Fox, L., Henrici, P., & Moler, C., Approximations and bounds for
eigenvalues of elliptic operators, SIAM Journal on Numerical Analysis,
4(1), 89–102 (1967).

Betcke, T., & Trefethen, L. N., Reviving the method of particular
solutions, SIAM review, 47(3), 469–491 (2005).
