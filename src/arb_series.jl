"""
    bessel_j(ν::arb, z::arb_series[, n = length(z)])

> Compute the Taylor series of the Bessel function.

It's computed using a recursive formula for the Taylor coefficients
of the Bessel functions and then composing it with `z`.
"""
function bessel_j(ν::arb, z::arb_series, n = length(z))
    res = arb_series(parent(z.poly)(), n)

    if n > 0
        x = z[0]
        a0 = bessel_j(ν, x)
        res[0] = a0
    end

    if n > 1
        a1 = 1//2*(bessel_j(ν - 1, x) - bessel_j(ν + 1, x))
        res[1] = a1
    end

    if n > 2
        a2 = 1//8*(bessel_j(ν - 2, x) + bessel_j(ν + 2, x) - 2a0)
        res[2] = a2
    end

    if n > 3
        a3 = 1//48*(bessel_j(ν - 3, x) - bessel_j(ν + 3, x) - 6a1)
        res[3] = a3
    end

    for i in 4:n-1
        k = i - 4
        ai = -(res[k]
               + 2x*res[k + 1]
               + (k^2 - ν^2 + x^2 + 4k + 4)*res[k + 2]
               + (2k^2 + 11k + 15)*x*res[k + 3]
               )/(x^2*(k^2 + 7k + 12))
        res[i] = ai
    end

    # Compose the Taylor series for the Bessel function with that of z
    z_tmp = arb_series(deepcopy(z.poly))
    z_tmp[0] = base_ring(parent(z.poly))(0)

    return Nemo.compose(res, z_tmp, n)
end

"""
    legendre_p(ν::arb, μ::arb, z::arb_series[, n = length(z)])

> Compute the Taylor series of the Legendre function.

It's computed using a recursive formula for the Taylor coefficients
of the Legendre functions and then composing it with `z`.
"""
function legendre_p(ν::arb, μ::arb, z::arb_series, n = length(z))
    # TODO: Precompute values used several times
    res = arb_series(parent(z.poly)(), n)

    if n > 0
        x = z[0]
        a0 = legendre_p_safe(ν, μ, x)
        res[0] = a0
    end

    if n > 1
        a1 = ((ν + 1)*x*a0 - (ν - μ + 1)*legendre_p_safe(ν + 1, μ, x))/(1 - x^2)
        res[1] = a1
    end

    if n > 2
        a2 = (2x*a1 - (ν*(ν + 1) - μ^2/(1 - x^2))*a0)/(2(1 - x^2))
        res[2] = a2
    end

    if n > 3
        a3 = (4(1 - x^2)^2*x*a2
              + 2(2μ^2 - (ν^2 + ν)*(1 - x^2))*x*a0
              + ((ν^2 + ν + 2)*x^2 + μ^2 - ν^2 - ν + 2)*(1 - x^2)*a1)/(6*(1 - x^2)^3)
        res[3] = a3
    end

    for i in 4:n-1
        k = i - 4
        ai = ((k + 1 + ν)*(ν - k)*res[k]
              - 4*x*(k^2 + (5k + 3 - ν^2 - ν)/2)*res[k + 1]
              + ((-6k^2 - 24k - 24 + ν^2 + ν)*x^2 + 2k^2 + 8k + 8 + μ^2 - ν^2 - ν)*res[k + 2]
              - 2(k + 3)*(2k + 5)*(x + 1)*x*(x - 1)*res[k + 3]
              )/((x - 1)^2*(x + 1)^2*(k + 4)*(k + 3))
        res[i] = ai
    end

    # Compose the Taylor series for the Legendre function with that of
    # z
    z_tmp = arb_series(deepcopy(z.poly))
    z_tmp[0] = base_ring(parent(z.poly))(0)

    return Nemo.compose(res, z_tmp, n)
end

function legendre_p_safe(ν::arb, μ::arb, z::arb_series, n = length(z))
    legendre_p(ν, μ, z, n)
end
