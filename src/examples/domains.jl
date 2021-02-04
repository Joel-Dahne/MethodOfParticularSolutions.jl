function example_domain_H(parent = RealField(precision(BigFloat)))
    angles =
        fmpq.([
            3 // 2,
            1 // 2,
            1 // 2,
            1 // 2,
            1 // 2,
            3 // 2,
            3 // 2,
            1 // 2,
            1 // 2,
            1 // 2,
            1 // 2,
            3 // 2,
        ])
    vertices = [
        parent.(v)
        for
        v in [
            (0, 0),
            (0, 1),
            (-1, 1),
            (-1, -2),
            (0, -2),
            (0, -1),
            (1, -1),
            (1, -2),
            (2, -2),
            (2, 1),
            (1, 1),
            (1, 0),
        ]
    ]

    domain = Polygon(angles, vertices, parent)

    u1 = StandaloneVertexEigenfunction(vertex(domain, 1), fmpq(1 // 2), fmpq(3 // 2))
    u2 = StandaloneVertexEigenfunction(vertex(domain, 6), fmpq(0), fmpq(3 // 2))
    u3 = StandaloneVertexEigenfunction(vertex(domain, 7), fmpq(-1 // 2), fmpq(3 // 2))
    u4 = StandaloneVertexEigenfunction(vertex(domain, 12), fmpq(1), fmpq(3 // 2))


    u = CombinedEigenfunction(
        domain,
        [u1, u2, u3, u4],
        us_to_boundary = [
            setdiff(boundaries(domain), (12, 1)),
            setdiff(boundaries(domain), (5, 6)),
            setdiff(boundaries(domain), (6, 7)),
            setdiff(boundaries(domain), (11, 12)),
        ],
    )

    return domain, u
end

function example_domain_ngon(
    n,
    parent = RealField(precision(BigFloat));
    lightning = false,
    linked = false,
)
    θ = fmpq((n - 2) // n)
    angles = fill(θ, n)
    vertices = [(cos(θ), sin(θ)) for θ in (2parent(π) / n) .* (0:n-1)]
    domain = Polygon(angles, vertices, parent)

    if !lightning
        if linked
            us = [LinkedEigenfunction(
                [
                    StandaloneVertexEigenfunction(
                        vertex(domain, i),
                        i * (1 - θ) + θ * 1 // 2,
                        θ,
                    ) for i in boundaries(domain)
                ],
                excluded_boundaries = [
                    BitSet([mod1(i - 1, n), i]) for i in boundaries(domain)
                ],
            )]
            us_to_boundary = fill(BitSet([1]), length(us))
        else
            us = [
                StandaloneVertexEigenfunction(
                    vertex(domain, i),
                    i * (1 - θ) + θ * 1 // 2,
                    θ,
                ) for i in boundaries(domain)
            ]
            us_to_boundary = [
                setdiff(boundaries(domain), (i, mod1(i - 1, length(boundaries(domain))))) for i in eachindex(us)
            ]
        end
    else
        if linked
            us = [
                LinkedEigenfunction(
                    [
                        StandaloneLightningEigenfunction(
                            vertex(domain, i),
                            i * (1 - θ) + θ * 1 // 2,
                            θ,
                        ) for i in boundaries(domain)
                    ],
                    ones(n),
                ),
                StandaloneInteriorEigenfunction(domain, stride = n),
            ]

            us_to_boundary = fill(BitSet([1]), length(us))
        else
            us = [
                [
                    StandaloneLightningEigenfunction(
                        vertex(domain, i),
                        i * (1 - θ) + θ * 1 // 2,
                        θ,
                    ) for i in boundaries(domain)
                ]
                StandaloneInteriorEigenfunction(domain)
            ]

            us_to_boundary = fill(boundaries(domain), length(us))
        end
    end

    u = CombinedEigenfunction(domain, us, us_to_boundary = us_to_boundary)

    return domain, u
end

function example_domain_triangle_in_triangle(parent = RealField(precision(BigFloat)))
    domain1 = Triangle(fmpq(1 // 3), fmpq(1 // 3), parent)
    domain2 = TransformedDomain(
        Triangle(fmpq(1 // 3), fmpq(1 // 3), parent),
        fmpq(0),
        parent(0.25),
        SVector(cospi(fmpq(1 // 6), parent), sinpi(fmpq(1 // 6), parent)) / 4,
    )

    domain = IntersectedDomain(domain1, domain2)

    u1 = StandaloneLightningEigenfunction(domain1, 1)

    u2 = LinkedEigenfunction([StandaloneLightningEigenfunction(domain1, i) for i = 2:3])

    u3 = StandaloneLightningEigenfunction(domain2, 1, l = parent(0.15), outside = true)

    u4 = LinkedEigenfunction([
        StandaloneLightningEigenfunction(domain2, i, l = parent(0.15), outside = true)
        for i = 2:3
    ])

    u5 = StandaloneInteriorEigenfunction(domain)


    us = [u1, u2, u3, u4, u5]
    #us = vcat([u1], u2.us, [u3], u4.us, [u5])
    orders = [1, 1, 1, 1, 1]
    us_to_boundary = fill(1:6, length(us))

    u = CombinedEigenfunction(domain, us, orders, us_to_boundary = us_to_boundary)

    return domain, u
end

function example_domain_ngon_in_ngon(
    n1,
    n2,
    parent = RealField(precision(BigFloat));
    linked = true,
    even = true,
    T = arb,
)
    θ1 = fmpq((n1 - 2) // n1)
    angles1 = fill(θ1, n1)
    vertices1 = [(cos(θ), sin(θ)) for θ in (2parent(π) / n1) .* (0:n1-1)]
    domain1 = Polygon(angles1, vertices1, parent)

    θ2 = fmpq((n2 - 2) // n2)
    angles2 = fill(θ2, n2)
    vertices2 = [(cos(θ), sin(θ)) for θ in (2parent(π) / n2) .* (0:n2-1)]
    domain2 = TransformedDomain(
        Polygon(angles2, vertices2, parent),
        fmpq(0),
        parent(0.5),
        SVector(parent(0), parent(0)),
    )

    domain = IntersectedDomain(domain1, domain2)

    if linked
        us = [
            LinkedEigenfunction([
                StandaloneLightningEigenfunction{T,arb}(domain1, i; even)
                for i in boundaries(domain1)
            ]),
            LinkedEigenfunction([
                StandaloneLightningEigenfunction{T,arb}(
                    domain2,
                    i,
                    outside = true,
                    l = parent(0.4);
                    even,
                ) for i in boundaries(domain2)
            ]),
            StandaloneInteriorEigenfunction(domain, stride = 4; even),
        ]

        orders = ifelse(even, [2, 2, 1], [3, 3, 2])
        even_boundaries = ifelse(even, [1, 5], [])

        u = CombinedEigenfunction(
            domain,
            us,
            orders,
            us_to_boundary = fill(BitSet([1, 5]), length(us));
            even_boundaries,
        )
    else
        us = vcat(
            [
                StandaloneLightningEigenfunction{T,arb}(domain1, i)
                for i in boundaries(domain1)
            ],
            [
                StandaloneLightningEigenfunction{T,arb}(
                    domain2,
                    i,
                    outside = true,
                    l = parent(0.4),
                ) for i in boundaries(domain2)
            ],
            [StandaloneInteriorEigenfunction(domain)],
        )

        u = CombinedEigenfunction(domain, us, [3, 3, 3, 3, 3, 3, 3, 3, 2])
    end

    return domain, u
end

"""
    example_domain_goal_v1(parent = RealField(precision(BigFloat)))

Return the domain given by a hexagon with 6 triangles cut out as well
as a corresponding eigenfunction.

Computing the first few eigenvalues in Matlab using a finite element
method gives (with `(N, d, h) = (22, 9, 6)`)
```
   31.5943
   75.8995
   75.9633
   75.9633
```
though it's likely not accurate to more than 1 or 2 decimals. If `(N,
d, h) = (27, 11, 6)` then they are
```
   31.0432
   63.2104
   63.7259
   63.7259
```
"""
function example_domain_goal_v1(
    N::Integer = 27,
    d::Integer = 11,
    h::Integer = 6,
    parent = RealField(precision(BigFloat));
    T = arb,
    lightning = false,
    inner_expansion = true,
    outer_expansion = false,
    symmetry_class = 1,
    even = true,
    reversed = true,
)
    ### Compute the domain ###
    # The main domain is a hexagon
    n = 6
    θ = fmpq((n - 2) // n)
    angles = fill(θ, n)
    vertices = [(cospi(θ, parent), sinpi(θ, parent)) for θ in fmpq(2 // n) .* (0:n-1)]
    exterior = Polygon(angles, vertices, parent)

    # The interior domains are triangles
    @assert mod(h, 3) == 0
    points = [
        SVector(parent(d // N), -parent(h // N) / sqrt(parent(3))),
        SVector(parent((d + h) // N), parent(0)),
        SVector(parent(d // N), parent(h // N) / sqrt(parent(3))),
    ]

    interior_angles = [1 // 3, 1 // 3, 1 // 3]

    interiors = [
        TransformedDomain(Polygon(interior_angles, points, parent), i // 3, 1, [0, 0])
        for i = 0:5
    ]

    domain = IntersectedDomain(exterior, interiors)

    ### Compute the eigenfunction ###
    # Expansions from the vertices of the hexagon
    if symmetry_class == 1
        if lightning
            outer_eigenfunction =
                i -> StandaloneLightningEigenfunction(
                    vertex(exterior, i),
                    fmpq(mod(Rational(1 - θ // 2 + (i - 1) * (1 - θ)), 2)),
                    θ;
                    even,
                )
            outer_excluded_boundary = i -> BitSet()
            order1 = ifelse(even, 2, 3)

        else
            outer_eigenfunction =
                i -> StandaloneVertexEigenfunction(
                    vertex(exterior, i),
                    i * (1 - θ) + θ * 1 // 2,
                    θ,
                    ifelse(even, 2, 1),
                )
            outer_excluded_boundary =
                i -> BitSet([mod1(i - 1, length(boundaries(exterior))), i])
            order1 = ifelse(even, 1, 2)
        end

        u1 = LinkedEigenfunction(
            [outer_eigenfunction(i) for i in boundaries(exterior)],
            excluded_boundaries = [
                outer_excluded_boundary(i) for i in boundaries(exterior)
            ],
        )

        # Expansions from the outer tip of the triangles
        u2 = LinkedEigenfunction([
            StandaloneLightningEigenfunction{T,fmpq}(
                vertex(d, 2),
                7 // 6 + d.rotation,
                2 - d.original.angles[2],
                l = parent(h // 2N),
                even = even,
            ) for d in interiors
        ])
        if even && reversed
            order2 = 3 * 2
        elseif even
            order2 = 2
        else
            order2 = 3
        end

        # Expansions from the inner tips of the triangles
        u3 = LinkedEigenfunction([
            StandaloneLightningEigenfunction{T,fmpq}(
                vertex(d, i),
                ifelse(i == 1, 1 // 2, -1 // 6) + d.rotation,
                2 - d.original.angles[i],
                l = parent(h // 2N),
                even = even && !reversed,
                reversed = reversed && i == 3,
            ) for d in interiors, i in [1, 3]
        ][:])
        if even && reversed
            order3 = 3 * 3
        elseif even
            order3 = 2
        else
            order3 = 3
        end

        # Expansion from the center
        u4 = StandaloneInteriorEigenfunction(domain, stride = 6; even)
        order4 = ifelse(even, 1, 2)

        # Interior expansions closer to outer boundaries
        u5 = LinkedEigenfunction([
            StandaloneInteriorEigenfunction(
                parent(3 // 4) .*
                SVector(cospi(fmpq(i // 3), parent), sinpi(fmpq(i // 3), parent)),
                fmpq(i // 3),
                stride = 6,
                even = true,
            ) for i = 0:5
        ])
        order5 = ifelse(even, 1, 2)

        us = AbstractPlanarEigenfunction[u1, u2, u3]
        orders = Int[order1, order2, order3]

        if inner_expansion
            push!(us, u4)
            push!(orders, order4)
        end
        if outer_expansion
            push!(us, u5)
            push!(orders, order5)
        end

        us_to_boundary = fill(BitSet([1, 7, 9]), length(us))
        even_boundaries = ifelse(even, Int[1, 9], Int[])
    elseif symmetry_class == 2
        ####################################################################
        is_even = i -> i == 1 || i == 4
        is_reversed = i -> i == 3 || i == 6

        if lightning
            outer_eigenfunction =
                i -> StandaloneLightningEigenfunction(
                    vertex(exterior, i),
                    fmpq(mod(Rational(1 - θ // 2 + (i - 1) * (1 - θ)), 2)),
                    θ;
                    even = is_even(i),
                    reversed = is_reversed(i),
                )
            outer_excluded_boundary = i -> BitSet()
            order11 = 2
            order12 = 3
        else
            outer_eigenfunction =
                i -> StandaloneVertexEigenfunction(
                    vertex(exterior, i),
                    i * (1 - θ) + θ * 1 // 2,
                    θ,
                    stride = ifelse(is_even(i), 2, 1),
                    reversed = is_reversed(i),
                )
            outer_excluded_boundary =
                i -> BitSet([mod1(i - 1, length(boundaries(exterior))), i])
            order11 = 1
            order12 = 1
        end

        u11 = LinkedEigenfunction(
            [outer_eigenfunction(i) for i in (1, 4)],
            [parent(1), parent(-1)],
            excluded_boundaries = [outer_excluded_boundary(i) for i in (1, 4)],
        )

        u12 = LinkedEigenfunction(
            [outer_eigenfunction(i) for i in (2, 3, 5, 6)],
            [parent(1), parent(-1), parent(-1), parent(1)],
            excluded_boundaries = [outer_excluded_boundary(i) for i in (2, 3, 5, 6)],
        )

        # Expansions from the outer tip of the triangles
        u21 = LinkedEigenfunction(
            [
                StandaloneLightningEigenfunction{T,fmpq}(
                    vertex(d, 2),
                    7 // 6 + d.rotation,
                    2 - d.original.angles[2],
                    l = parent(h // 2N),
                    even = is_even(i),
                ) for (i, d) in ((1, interiors[1]), (4, interiors[4]))
            ],
            [parent(1), parent(-1)],
        )
        order21 = 2

        u22 = LinkedEigenfunction(
            [
                StandaloneLightningEigenfunction{T,fmpq}(
                    vertex(d, 2),
                    7 // 6 + d.rotation,
                    2 - d.original.angles[2],
                    l = parent(h // 2N),
                    even = is_even(i),
                    reversed = is_reversed(i),
                )
                for
                (i, d) in (
                    (2, interiors[2]),
                    (3, interiors[3]),
                    (5, interiors[5]),
                    (6, interiors[6]),
                )
            ],
            [parent(1), parent(-1), parent(-1), parent(1)],
        )
        order22 = 3

        # Expansions from the inner tips of the triangles
        u31 = LinkedEigenfunction(
            [
                StandaloneLightningEigenfunction{T,fmpq}(
                    vertex(d, i),
                    ifelse(i == 1, 1 // 2, -1 // 6) + d.rotation,
                    2 - d.original.angles[i],
                    l = parent(h // 2N),
                    reversed = i == 3,
                ) for (j, d) in ((1, interiors[1]), (4, interiors[4])), i in [1, 3]
            ][:],
            [parent(1), parent(-1), parent(1), parent(-1)],
        )
        order31 = 3

        u32 = LinkedEigenfunction(
            [
                StandaloneLightningEigenfunction{T,fmpq}(
                    vertex(d, i),
                    ifelse(i == 1, 1 // 2, -1 // 6) + d.rotation,
                    2 - d.original.angles[i],
                    l = parent(h // 2N),
                    reversed = i == 3,
                )
                for
                (j, d, i) in (
                    (2, interiors[2], 1),
                    (3, interiors[3], 3),
                    (5, interiors[5], 1),
                    (6, interiors[6], 3),
                )
            ],
            [parent(1), parent(-1), parent(-1), parent(1)],
        )
        order32 = 3

        u33 = LinkedEigenfunction(
            [
                StandaloneLightningEigenfunction{T,fmpq}(
                    vertex(d, i),
                    ifelse(i == 1, 1 // 2, -1 // 6) + d.rotation,
                    2 - d.original.angles[i],
                    l = parent(h // 2N),
                    reversed = i == 3,
                )
                for
                (j, d, i) in (
                    (2, interiors[2], 3),
                    (3, interiors[3], 1),
                    (5, interiors[5], 3),
                    (6, interiors[6], 1),
                )
            ],
            [parent(1), parent(-1), parent(-1), parent(1)],
        )
        order33 = 3

        # Expansion from the center
        u4 = StandaloneInteriorEigenfunction(domain, stride = 2, offset = 1; even)
        order4 = 1

        us = AbstractPlanarEigenfunction[u11, u12, u21, u22, u31, u32, u33, u4]
        orders =
            Int[order11, order12, 3order21, 3order22, 3order31, 3order32, 3order33, order4]

        us_to_boundary = fill(BitSet([1, 2, 7, 9, 10, 11, 12]), length(us))
        even_boundaries = Int[2, 9]
    else
        throw(ArgumentError("symmetry_class should be 1 or 2, got $symmetry_class"))
    end

    u = CombinedEigenfunction(domain, us, orders; us_to_boundary, even_boundaries)

    return domain, u
end

"""
    example_domain_goal_v2(parent = RealField(precision(BigFloat)))

Return the domain given by a hexagon with 6 polygon cut out as well
as a corresponding eigenfunction.

With `a = 0.18`, `b = 0.1` and `c = 0.2` computing the first few
eigenvalues in Matlab using a finite element method gives
```
   26.3259
   56.8914
   56.9041
   56.9041
```
If `rotated` is true they are
```
   26.3348
   59.1647
   59.3959
   59.3959
```
"""
function example_domain_goal_v2(
    parent = RealField(precision(BigFloat));
    T = arb,
    even = true,
    reversed = true,
    inner_expansion = true,
    outer_expansion = false,
    a = 0.18,
    b = 0.1,
    c = 0.2,
    rotated = true,
)
    # The main domain is a hexagon
    n = 6
    θ = fmpq((n - 2) // n)
    angles = fill(θ, n)
    vertices = [(cospi(θ, parent), sinpi(θ, parent)) for θ in fmpq(2 // n) .* (0:n-1)]
    exterior = Polygon(angles, vertices, parent)

    # Points for the interior domains
    points = [
        SVector(parent(1), parent(0)),
        SVector(parent(1), parent(0)),
        SVector(cospi(fmpq(1 // 3), parent), sinpi(fmpq(1 // 3), parent)),
        SVector(cospi(fmpq(1 // 3), parent), sinpi(fmpq(1 // 3), parent)),
    ]

    points[1] *= a + b
    points[2] *= a
    points[3] *= a
    points[4] *= a + b

    points[2] +=
        2 * sinpi(fmpq(1 // 3), parent) * b .*
        SVector(cospi(fmpq(1 // 6), parent), sinpi(fmpq(1 // 6), parent))
    points[3] +=
        2 * sinpi(fmpq(1 // 3), parent) * b .*
        SVector(cospi(fmpq(1 // 6), parent), sinpi(fmpq(1 // 6), parent))

    interior_angles = [1 // 3, 2 // 3, 2 // 3, 1 // 3]

    ϕ = ifelse(rotated, -1 // 6, 0 // 6)

    interiors = [
        TransformedDomain(
            Polygon(interior_angles, points, parent),
            fmpq(i // 3 + ϕ),
            parent(1),
            c .* SVector(
                cospi(fmpq(i // 3 + 1 // 6 + ϕ), parent),
                sinpi(fmpq(i // 3 + 1 // 6 + ϕ), parent),
            ),
        ) for i = 0:5
    ]

    domain = IntersectedDomain(exterior, interiors)

    # Expansions from the vertices of the hexagon
    u1 = LinkedEigenfunction([
        StandaloneLightningEigenfunction{arb,fmpq}(
            vertex(exterior, i),
            fmpq(mod(Rational(1 - θ // 2 + (i - 1) * (1 - θ)), 2)),
            θ;
            even,
        ) for i in boundaries(exterior)
    ])

    # Expansions from the outer part of the interior polygons
    u2 = LinkedEigenfunction([
        StandaloneLightningEigenfunction{T,fmpq}(
            vertex(d, i),
            ifelse(i == 2, fmpq(4 // 3), fmpq(5 // 3)) + fmpq(d.rotation),
            2 - d.original.angles[i],
            l = parent(0.08),
            even = even && !reversed,
            reversed = reversed && i == 3,
        ) for d in interiors, i in [2, 3]
    ][:])

    # Expansions from the inner part of the interior polygons
    u3 = LinkedEigenfunction([
        StandaloneLightningEigenfunction{T,fmpq}(
            vertex(d, i),
            ifelse(i == 1, fmpq(2 // 3), fmpq(0)) + fmpq(d.rotation),
            2 - d.original.angles[i],
            l = parent(0.08),
            even = even && !reversed,
            reversed = reversed && i == 4,
        ) for d in interiors, i in [1, 4]
    ][:])

    # Expansion from the center
    u4 = StandaloneInteriorEigenfunction(domain, stride = 6; even)

    # Interior expansions closer to outer boundaries
    u5 = LinkedEigenfunction([
        StandaloneInteriorEigenfunction(
            3 // 4 .* SVector(cospi(fmpq(i // 3), parent), sinpi(fmpq(i // 3), parent)),
            fmpq(i // 3),
            stride = 6,
            even = true,
        ) for i = 0:5
    ])

    us = AbstractPlanarEigenfunction[u1, u2, u3]

    if even && reversed
        orders = [2, 3, 3]
    elseif even
        orders = [2, 2, 2]
    else
        orders = [3, 3, 3]
    end

    if inner_expansion
        push!(us, u4)
        push!(orders, ifelse(even, 1, 2))
    end
    if outer_expansion
        push!(us, u5)
        push!(orders, ifelse(even, 1, 2))
    end

    us_to_boundary = fill(BitSet([1, 7, 8, 10]), length(us))
    even_boundaries = ifelse(even, Int[1, 8, 10], Int[])

    u = CombinedEigenfunction(domain, us, orders; us_to_boundary, even_boundaries)

    return domain, u
end
