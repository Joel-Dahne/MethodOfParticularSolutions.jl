# MethodOfParticularSolutions.jl

This is a package for computing eigenvalues and eigenfunctions of the
Laplacian on planar or spherical domains using the method of
particular solutions in Julia.

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

## References

Fox, L., P. Henrici, and C. Moler. "Approximations and bounds for eigenvalues of elliptic operators." SIAM Journal on Numerical Analysis 4.1 (1967): 89-102.

Betcke, Timo, and Lloyd N. Trefethen. "Reviving the method of particular solutions." SIAM review 47.3 (2005): 469-491.
