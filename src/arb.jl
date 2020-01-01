"""
    atan(x::arb, y::arb)
> Return atan(x, y) = arg(x + yi).
"""
function atan(x::arb, y::arb)
    z = parent(x)()
    ccall((:arb_atan2, Nemo.libarb), Nothing,
          (Ref{arb}, Ref{arb}, Ref{arb}, Int), z, x, y, parent(x).prec)
    return z
end

"""
    legendre_p(ν::arb, μ::arb, z::arb)
> Returns the associated Legendre function of the first kind evaluated
  for degree ν, order μ, and argument z. When μ is zero, this reduces
  to the Legendre polynomial P_ν(z).
"""
function legendre_p(ν::arb, μ::arb, z::arb)
    res = parent(z)()
    ccall((:arb_hypgeom_legendre_p, Nemo.libarb), Nothing,
          (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong),
          res, ν, μ, z, 0, parent(z).prec)
    return res
end

function legendre_p_safe(ν::arb, μ::arb, z::arb)
    res = parent(z)(NaN)
    prec = parent(z).prec
    while !isfinite(res) && prec <= 8*parent(z).prec
        ccall((:arb_hypgeom_legendre_p, Nemo.libarb), Nothing,
              (Ref{arb}, Ref{arb}, Ref{arb}, Ref{arb}, Cint, Clong),
              res, ν, μ, z, 0, prec)
        prec *= 2
    end

    return res
end

"""
    getinterval(x::arb)
> Return an interval [a,b] containing the ball x.
"""
function getinterval(x::arb)
    getinterval(arb, x)
end

"""
    getinterval(::Type{arb}, x::arb)
> Return an interval [a,b] containing the ball x.
"""
function getinterval(::Type{arb}, x::arb)
    a, b = x.parent(), x.parent()
    a_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), a)
    b_mid = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), b)
    ccall((:arb_get_interval_arf, Nemo.libarb), Cvoid,
          (Ptr{Nemo.arf_struct}, Ptr{Nemo.arf_struct}, Ref{arb}, Clong),
          a_mid, b_mid, x, x.parent.prec)

    (a, b)
end

"""
    getinterval(::Type{BigFloat}, x::arb)
> Return an interval [a,b] containing the ball x.
"""
function getinterval(::Type{BigFloat}, x::arb)
    a, b = BigFloat(), BigFloat()
    ccall((:arb_get_interval_mpfr, Nemo.libarb), Cvoid,
          (Ref{BigFloat}, Ref{BigFloat}, Ref{arb}),
          a, b, x)

    (a, b)
end

"""
    convert(::Type{BigFloat}, x::arb)
> Return the midpoint of x as a BigFloat rounded down to the current
  precision of BigFloat.
"""
function BigFloat(x::arb)
    GC.@preserve x begin
        t = ccall((:arb_mid_ptr, Nemo.libarb), Ptr{Nemo.arf_struct}, (Ref{arb}, ), x)
        # 4 == round to nearest
        m = BigFloat()
        ccall((:arf_get_mpfr, Nemo.libarb), Float64,
              (Ref{BigFloat}, Ptr{Nemo.arf_struct}, Base.MPFR.MPFRRoundingMode),
              m, t, Base.MPFR.MPFRRoundNearest)
    end
    return m
end

"""
    convert(::Type{BigFloat}, x::arb)
> Return the midpoint of x as a BigFloat rounded down to the current
  precision of BigFloat.
"""
function Base.convert(::Type{BigFloat}, x::arb)
    return BigFloat(x)
end
