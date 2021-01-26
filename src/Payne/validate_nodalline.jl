"""
    validate_nodalline(domain, u, λ, distance)

Return true if `u` is strictly positive/negative on a loop and
negative/positive on at least one point inside.

The loop is given by line between `(distance, 0)` and `(distance,
tan(π/6)*distance)` which is extended by symmetry to loop around the
midpoint of the domain.

FIXME: Doesn't currently take the error of `u` into account. We have
to find the L-inf error.
"""
function validate_nodalline(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    λ::arb,
    distance::arb,
)
    # Check that the line indeed is inside the interior domains
    @assert distance < vertex(domain.interiors[1], 1)[1]

    # First step is to check the sign of u inside
    inside = u(SVector(distance/2, zero(distance)), λ)

    if isnegative(inside)
        sign_inside = -1
    elseif ispositive(inside)
        sign_inside = 1
    else
        @error "evaluation on the interior point contains zero"
        @show inside
        return false
    end

    ### Second step is to prove that u has the opposite sign on a loop
    ### which encloses the middle but is inside all the interior
    ### domains. By symmetry arguments we can reduce this to checking
    ### u on the line between (distance, 0) and (distance,
    ### tan(π/6)*distance).

    # Parameterization of the line between (distance, 0) and
    # (distance, tan(π/6)*distance)
    y_max = distance*tanpi(domain.parent(1//6))
    parameterization(t::arb) =  SVector(distance, t*y_max)
    parameterization(t::arb_series) =  SVector(distance + 0*t, t*y_max)

    # Function for computing u along this line. Multiply u by the
    # opposite sign of the inverse so that it's supposed to be
    # negative.
    @show sign_inside
    if sign_inside < 0
        f = t -> -u(parameterization(t), λ)
    else
        f = t -> u(parameterization(t), λ)
    end

    @show f(zero(λ))

    # Prove that u multiplied by the opposite sign of the inverse is
    # strictly negative.
    strictly_positive, enclosure = bounded_by(
        f,
        zero(λ),
        one(λ),
        zero(λ),
        use_taylor = true,
        n = length(coefficients(u)),
        show_trace = true,
        return_enclosure = true,
    )

    @show enclosure

    if ismissing(strictly_positive)
        @error "could not prove that u was strictly positive on the line"
        return false
    elseif !strictly_positive
        @error "u is not strictly positive on the line"
        return false
    else
        @info "Managed to prove that u is strictly positive on the line"
    end

    return true
end
