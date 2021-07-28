function eigenfunction_plotdata(u::StandaloneVertexEigenfunction)
    v = convert(SVector{2,Float64}, u.vertex)
    if has_rational_angles(u)
        # We need Rational here since fmpq doesn't convert directly to
        # Float64
        orientation = π * convert(Float64, Rational(u.orientation))
        θ = π * convert(Float64, Rational(u.θ))
    else
        orientation = convert(Float64, u.orientation)
        θ = convert(Float64, u.θ)
    end

    vertex = ([v[1]], [v[2]])
    edges = begin
        v1 = v + 0.1 * [cos(orientation), sin(orientation)]
        v2 = v + 0.1 * [cos((orientation + θ)), sin((orientation + θ))]
        ([v1[1], v[1], v2[1]], [v1[2], v[2], v2[2]])
    end
    arc = begin
        res = [
            v + 0.1 * [cos(orientation + t * θ), sin(orientation + t * θ)]
            for t in range(0, 1, length = 20)
        ]
        (getindex.(res, 1), getindex.(res, 2))
    end

    return vertex, edges, arc
end

function eigenfunction_plotdata(u::StandaloneInteriorEigenfunction)
    v = convert(SVector{2,Float64}, u.vertex)
    if has_rational_angles(u)
        # We need Rational here since fmpq doesn't convert directly to
        # Float64
        orientation = π * convert(Float64, Rational(u.orientation))
    else
        orientation = convert(Float64, u.orientation)
    end

    vertex = ([v[1]], [v[2]])
    edge = begin
        w = v + 0.1 * [cos(orientation), sin(orientation)]
        ([v[1], w[1]], [v[2], w[2]])
    end
    arc = begin
        res = [v + 0.1 * [cos(θ), sin(θ)] for θ in range(0, 2π, length = 50)]
        (getindex.(res, 1), getindex.(res, 2))
    end

    return vertex, edge, arc
end

function eigenfunction_plotdata(u::StandaloneLightningEigenfunction)
    v = convert(SVector{2,Float64}, u.vertex)
    if has_rational_angles(u)
        # We need Rational here since fmpq doesn't convert directly to
        # Float64
        orientation = π * convert(Float64, Rational(u.orientation))
        θ = π * convert(Float64, Rational(u.θ))
    else
        orientation = convert(Float64, u.orientation)
        θ = convert(Float64, u.θ)
    end

    vertex = ([v[1]], [v[2]])
    edges = begin
        v1 = v + 0.1 * [cos(orientation), sin(orientation)]
        v2 = v + 0.1 * [cos((orientation + θ)), sin((orientation + θ))]
        ([v1[1], v[1], v2[1]], [v1[2], v[2], v2[2]])
    end
    arc = begin
        res = [
            v + 0.1 * [cos(orientation + t * θ), sin(orientation + t * θ)]
            for t in range(0, 1, length = 20)
        ]
        (getindex.(res, 1), getindex.(res, 2))
    end
    charges = begin
        n = div(length(coefficients(u)) - 1, 3) + 1
        cs = [charge(u, i, n, true) for i = 1:n]
        (Float64.(getindex.(cs, 1)), Float64.(getindex.(cs, 2)))
    end

    return vertex, edges, arc, charges
end

@recipe function f(u::StandaloneVertexEigenfunction)
    seriestype := [:path :path :scatter]
    linewidth := [2 3 3]
    markersize := 5
    label --> ["" "u" ""]
    seriescolor := [:blue :blue :red]
    linestyle := [:dot :solid :solid]

    vertex, edges, arc = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [edges[1]; fill(NaN, length(arc[1]) - length(edges[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
    )
    y = hcat(
        arc[2],
        [edges[2]; fill(NaN, length(arc[2]) - length(edges[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
    )

    return (x, y)
end

@recipe function f(u::StandaloneInteriorEigenfunction)
    seriestype := [:path :path :scatter]
    linewidth := [2 2 3]
    markersize := 5
    label --> ["" "" "u"]
    seriescolor := [:blue :blue :red]
    linestyle := [:dot :solid :solid]
    vertex, edge, arc = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [edge[1]; fill(NaN, length(arc[1]) - length(edge[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
    )
    y = hcat(
        arc[2],
        [edge[2]; fill(NaN, length(arc[2]) - length(edge[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
    )
    return (x, y)
end

@recipe function f(u::StandaloneLightningEigenfunction)
    seriestype --> [:path :path :scatter :scatter]
    linewidth := [2 3 3 2]
    markersize := [5 5 5 2]
    label --> ["" "" "u" ""]
    seriescolor := [:blue :blue :red :green]
    linestyle := [:dot :solid :solid :solid]

    vertex, edges, arc, charges = eigenfunction_plotdata(u)

    x = hcat(
        arc[1],
        [edges[1]; fill(NaN, length(arc[1]) - length(edges[1]))],
        [vertex[1]; fill(NaN, length(arc[1]) - length(vertex[1]))],
        [charges[1]; fill(NaN, length(arc[1]) - length(charges[1]))],
    )
    y = hcat(
        arc[2],
        [edges[2]; fill(NaN, length(arc[2]) - length(edges[2]))],
        [vertex[2]; fill(NaN, length(arc[2]) - length(vertex[2]))],
        [charges[2]; fill(NaN, length(arc[2]) - length(charges[2]))],
    )

    return (x, y)
end

@userplot EigenfunctionPlot

@recipe function f(h::EigenfunctionPlot)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: CombinedEigenfunction)
        error(
            "EigenfunctionPlots should be given a CombinedEigenfunction. Got: $(typeof(h.args))",
        )
    end
    u = h.args[1]

    for v in u.us
        @series begin
            label := ""
            v
        end
    end

    return nothing
end

@userplot EigenfunctionHeatmap

@recipe function f(
    h::EigenfunctionHeatmap;
    absolute_value = false,
    twosided = !absolute_value,
    include_exterior = false,
)
    if !(length(h.args) == 3 || length(h.args) == 5) ||
       !(typeof(h.args[1]) <: AbstractDomain) ||
       !(typeof(h.args[2]) <: AbstractEigenfunction) ||
       !(typeof(h.args[3]) <: Union{Real,arb})
        throw(
            ArgumentError(
                "EigenfunctionHeatmap should be given a domain, an eigenfunction, an eigenvalue " *
                "and x, y. Got: $(typeof(h.args))",
            ),
        )
    end
    if length(h.args) == 5
        domain, u, λ, xs, ys = h.args[1:5]
    elseif length(h.args) == 3
        domain, u, λ = h.args[1:3]
        xs = 50
        ys = 50
    end

    vs = vertices(domain)
    if typeof(xs) <: Integer
        xs = range(extrema(Float64.(collect(getindex.(vs, 1))))..., length = xs)
    end
    if typeof(ys) <: Integer
        ys = range(extrema(Float64.(collect(getindex.(vs, 2))))..., length = ys)
    end

    pts = SVector.(xs', ys)
    res = similar(pts, Float64)
    Threads.@threads for i in eachindex(pts)
        if include_exterior || pts[i] ∈ domain
            res[i] = ifelse(absolute_value, abs, identity)(u(pts[i], λ))
        else
            res[i] = 0
        end
    end

    xlims := extrema(xs)
    ylims := extrema(ys)
    label := ""
    aspect_ratio --> :equal

    if twosided
        # Make the color bar symmetric around 0
        m = maximum(abs, res)
        clims --> (-m, m)
    end

    @series begin
        seriestype := :heatmap
        if twosided
            seriescolor --> :balance
        else
            seriescolor --> :matter
        end
        xs, ys, res
    end

    @series begin
        domain, 100, 0
    end

    return nothing
end
