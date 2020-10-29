@recipe function f(domain::AbstractPlanarDomain, n::Integer = 100)
    seriestype := :scatter
    markersize --> 2
    label --> ["Boundary" "Interior"]

    # TODO: Make it so that you can have different number of points
    # for the boundary and interior
    boundary = boundary_points(domain, n)[1]
    interior = interior_points(domain, n)

    return (
        Float64[getindex.(boundary, 1) getindex.(interior, 1)],
        Float64[getindex.(boundary, 2) getindex.(interior, 2)],
    )
end
