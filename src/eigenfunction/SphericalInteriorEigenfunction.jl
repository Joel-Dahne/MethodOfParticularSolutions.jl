function SphericalInteriorEigenfunction(domain::SphericalTriangle,
                                        θ::arb,
                                        ϕ::arb)
    SphericalInteriorEigenfunction(domain, θ, ϕ, arb[])
end

"""
    coordinate_transformation(u::SphericalInteriorEigenfunction)
> Return a coordinate transformation T which switches from spherical
  coordinates (θ, ϕ) to (θ', ϕ'), in the new coordinate system the
  point from which the eigenfunction originates is put at the north
  pole.

  This is computed by switching to cartesian coordinates, performing
  two rotations and then switching back to spherical coordinates.
"""
function coordinate_transformation(u::SphericalInteriorEigenfunction)
    # TODO: The performance could most likely be improved

    # Rotate by u.ϕ along the along the z-axis so that the ϕ value of
    # the point becomes zero, the rotate by u.θ along the y-axis so
    # that the point ends up on the north pole.
    L = LinearMap(RotYZ(-u.θ, -u.ϕ))
    T = (θ, ϕ) -> begin
        xyz = cartesian(θ, ϕ)
        xyz = L(xyz)
        (θ, ϕ) = spherical(xyz)
        return θ, ϕ
    end

    return T
end

"""
    u(θ::arb, ϕ::arb, λ::arb, k::Integer; notransform::Bool = false)
> Evaluate the k-th basis function for the eigenfunction with the
  given λ on the point given by (θ, ϕ). If notransform is true then
  do not perform a coordinate transform on θ and ϕ first, this assumes
  that they already given in the coordinate system used by u. See
  coordinate_transform for details about the coordinate transform
  used.
"""
function (u::SphericalInteriorEigenfunction)(θ::arb,
                                             ϕ::arb,
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false)
    ν = u.domain.parent(-0.5) + sqrt(u.domain.parent(0.25) + λ)
    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end
    if k == 1
        μ = u.domain.parent(0)
        return legendre_p_safe(ν, μ, cos(θ))
    elseif k % 2 == 0
        μ = u.domain.parent(div(k, 2))
        return legendre_p_safe(ν, μ, cos(θ))*cos(μ*ϕ)
    else
        μ = u.domain.parent(div(k, 2))
        return legendre_p_safe(ν, μ, cos(θ))*sin(μ*ϕ)
    end
end

function (u::SphericalInteriorEigenfunction)((θ, ϕ),
                                             λ::arb,
                                             k::Integer;
                                             notransform::Bool = false)
    u(θ, ϕ, λ, k, notransform = notransform)
end

"""
    u(θ::arb, ϕ::arb, λ::arb; notransform::Bool = false)
> Evaluate the eigenfunction with the given λ on the point given by
  (θ, ϕ). If notransform is true then do not perform a coordinate
  transform on θ and ϕ first, this assumes that they already given in
  the coordinate system used by u. See coordinate_transform for
  details about the coordinate transform used.
"""
function (u::SphericalInteriorEigenfunction)(θ::arb,
                                             ϕ::arb,
                                             λ::arb;
                                             notransform::Bool = false)
    res = θ.parent(0)

    if !notransform
        θ, ϕ = coordinate_transformation(u, θ, ϕ)
    end

    for k in 1:length(u.coefficients)
        res += u.coefficients[k]*u(θ, ϕ, λ, k, notransform = true)
    end

    res
end

function (u::SphericalInteriorEigenfunction)((θ, ϕ),
                                             λ::arb;
                                             notransform::Bool = false)
    u(θ, ϕ, λ)
end
