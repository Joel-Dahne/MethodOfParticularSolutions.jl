"""
    enclose_eigenvalue_radius(λ, domain, defect, norm)

Compute the radius `r` of a ball such that the interval `[λ - r, λ +
r]` contains an eigenvalue of the given domain.

Given `ϵ = sqrt(area(domain)) * defect / norm` there is an eigenvalue
`λₖ` of the domain for some `k` satisfying
```
abs(λ - λₖ) / λₖ <= ϵ
```
which can be rewritten as
```
λ / (1 + ϵ) <= λₖ <= λ / (1 - ϵ)
```
Taking
```
r = max(λ - λ / (1 + ϵ), λ / (1 - ϵ) - λ)
```
we are therefore guaranteed to have `λₖ ∈ [λ - r, λ + r]`. Further we
can notice that the maximum is always given by `λ / (1 - ϵ) - λ`.

The value of `defect` and `norm` should be such as computed by
[`defect`](@ref) and [`norm`](@ref). **TODO:** Implement these
methods.

For rigorous results the parameters should be given as types
supporting rigorous computations. If rigorous results are not required
then [`defect_estimate`](@ref) and [`norm_estimate`](@ref) can also be
used for computing estimates of `defect` and `norm`.

For references about the enclosures see

> Moler, C. B. and Payne, L. E., Bounds for eigenvalues and eigenvectors of symmetric operators
> _SIAM Journal on Numerical Analysis, 5(1), 64–70 (1968)_
> https://doi.org/10.1137/0705004

> Fox, L., Henrici, P. and Moler, C., Approximations and bounds for eigenvalues of elliptic operators
> _SIAM Journal on Numerical Analysis, 4(1), 89–102 (1967)_
> https://doi.org/10.1137/0704008
"""
function enclose_eigenvalue_radius(λ, domain, defect, norm)
    ϵ = sqrt(area(domain)) * defect / norm

    r = λ / (1 - ϵ) - λ

    return r
end
