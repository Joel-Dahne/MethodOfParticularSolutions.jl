@recipe function f(domain::AbstractPlanarDomain, nb::Integer = 100, ni::Integer = 100)
    seriestype := :scatter
    markersize --> 2
    label --> ["Boundary" "Interior"]

    boundary = boundary_points(domain, nb, distribution = :linear)[1]
    interior = interior_points(domain, ni)

    x = Float64.(hcat(
        [getindex.(boundary, 1); fill(NaN, max(0, ni - nb))],
        [getindex.(interior, 1); fill(NaN, max(0, nb - ni))]
    ))
    y = Float64.(hcat(
        [getindex.(boundary, 2); fill(NaN, max(0, ni - nb))],
        [getindex.(interior, 2); fill(NaN, max(0, nb - ni))]
    ))

    return x, y
end
