"""
    continuedfraction(x::arb)
Compute a continued fraction representation of `x`.
"""
function continuedfraction(x::arb)
    i = floor(x)

    if isint(i)
        return [i; continuedfraction(inv(x - i))]
    else
        return [i]
    end
end

"""
    convergents(expansion::Vector{T}) where {T <: Integer}
Compute the convergents given a continued faction expansion `as`.
"""
function convergents(as::Vector{T}) where {T<:Integer}
    hs = zeros(T, length(as))
    ks = zeros(T, length(as))

    hs[1] = as[1]
    ks[1] = one(T)

    hs[2] = as[1] * as[2] + 1
    ks[2] = as[2]

    for n = 3:length(as)
        hs[n] = as[n] * hs[n-1] + hs[n-2]
        ks[n] = as[n] * ks[n-1] + ks[n-2]
    end

    return hs .// ks
end

"""
    convergents(expansion::Vector{arb})
Compute the convergents given a continued faction expansion `as`. The
computations are done using ball arithmetic and it returns a list of
pairs corresponding to the nominator and denominator of the
convergents.
"""
function convergents(as::Vector{arb})
    hs = similar(as)
    ks = similar(as)

    hs[1] = as[1]
    ks[1] = one(as[1])

    hs[2] = as[1] * as[2] + 1
    ks[2] = as[2]

    for n = 3:length(as)
        hs[n] = as[n] * hs[n-1] + hs[n-2]
        ks[n] = as[n] * ks[n-1] + ks[n-2]
    end

    return collect(zip(hs, ks))
end
