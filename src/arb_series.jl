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

    if !isfinite(a0)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
    end

    if n > 1
        a1 = 1//2*(bessel_j(ν - 1, x) - bessel_j(ν + 1, x))
        res[1] = a1
    end

    if !isfinite(a0)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
    end

    if n > 2
        a2 = 1//8*(bessel_j(ν - 2, x) + bessel_j(ν + 2, x) - 2a0)
        res[2] = a2
    end

    if !isfinite(a0)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
    end

    if n > 3
        a3 = 1//48*(bessel_j(ν - 3, x) - bessel_j(ν + 3, x) - 6a1)
        res[3] = a3
    end

    if !isfinite(a0)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
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
    legendre_p(ν::arb_series, μ::arb, z::arb[, n = length(ν)])

> Compute the Taylor series of the Legendre function with respect to
  the parameter ν.
"""
function legendre_p(ν::arb_series, μ::arb, z::arb, n = length(ν))
    CC = ComplexField(z.parent.prec)
    PP = AcbPolyRing(CC, :x)
    νν = PP(ν.poly)
    μμ = PP(μ)
    zz = PP(z)

    a = -νν
    b = νν + 1
    c = 1 - μμ

    res = PP()

    ccall((:acb_hypgeom_2f1_series_direct, :libarb), Cvoid, (Ref{acb_poly}, Ref{acb_poly},
                                                             Ref{acb_poly}, Ref{acb_poly},
                                                             Ref{acb_poly}, Cint, Clong, Clong),
          res, a, b, c, 0.5*(1 - zz), 1, 2, CC.prec)

    realres = arb_series(ν.poly.parent([real(coeff(res, i)) for i in 0:degree(res)]), n)
    ((1 + z)/(1 - z))^(μ/2)*realres
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

    if !isfinite(a0)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
    end

    x2 = x^2
    onemx2 = 1 - x2

    if n > 1
        a1 = ((ν + 1)*x*a0 - (ν - μ + 1)*legendre_p_safe(ν + 1, μ, x))/(onemx2)
        res[1] = a1
    end

    if !isfinite(a1)
        return arb_series(parent(z.poly)(fill(NaN, n)), n)
    end

    ν2 = ν^2
    μ2 = μ^2

    if n > 2
        a2 = (2x*a1 - (ν2 + ν - μ2/(onemx2))*a0)/(2(onemx2))
        res[2] = a2
    end

    if n > 3
        a3 = (4(onemx2)^2*x*a2
              + 2(2μ2 - (ν2 + ν)*(onemx2))*x*a0
              + ((ν2 + ν + 2)*x2 + μ2 - ν2 - ν + 2)*(onemx2)*a1)/(6*(onemx2)^3)
        res[3] = a3
    end

    onemx22 = onemx2^2
    xonemx2 = x*onemx2
    ν2pν = ν2 + ν
    μ2mν2pν = μ2 - ν2pν
    fourx = 4x
    for i in 4:n-1
        k = i - 4
        ai = ((k + 1 + ν)*(ν - k)*res[k]
              - (0.5*(2k^2 + (5k + 3) - ν2pν))*fourx*res[k + 1]
              + (((-6k^2 - 24k - 24) + ν2pν)*x2 + (2k^2 + 8k + 8) + μ2mν2pν)*res[k + 2]
              + (2(k + 3)*(2k + 5))*xonemx2*res[k + 3]
              )/(onemx22*((k + 4)*(k + 3)))
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
