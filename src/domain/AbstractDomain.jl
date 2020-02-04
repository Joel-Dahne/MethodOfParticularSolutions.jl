"""
    boundary_parameterization(t,
                              domain::AbstractDomain,
                              i::Integer)
> Compute a parameterization of boundary number `i`.

  The parameterization goes from t = 0 to t = 1. The enumeration of
  the boundary depends on the type of domain used.
"""
function boundary_parameterization(t,
                                   domain::AbstractDomain,
                                   i::Integer)
    throw(ErrorException("boundary_parameterization not implemented for domain of type $(typeof(domain))"))
end

"""
    boundary_points(domain::AbstractDomain,
                    u::AbstractEigenfunction,
                    n::Integer)
> Return n points from boundaries of the domain on which the
  eigenfunction is not identically zero.
"""
function boundary_points(domain::AbstractDomain,
                         u::AbstractEigenfunction,
                         n::Integer)
    throw(ErrorException("boundary_points not implemented for domain of type $(typeof(domain))"))
end

"""
    interior_points(domain::AbstractDomain,
                    n::Integer;
                    rng = MersenneTwister(42))
> Return n points taken uniformly at random from the interior of the
  domain.
"""
function interior_points(domain::AbstractDomain,
                         n::Integer;
                         rng = MersenneTwister(42))
    throw(ErrorException("interior_points not implemented for domain of type $(typeof(domain))"))
end
