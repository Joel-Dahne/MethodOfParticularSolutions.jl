function example_domain_H(parent = RealField(precision(BigFloat)))
    angles = fmpq.([3//2, 1//2, 1//2, 1//2, 1//2, 3//2, 3//2, 1//2, 1//2, 1//2, 1//2, 3//2])
    vertices = [parent.(v) for v in [
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
    ]]

    domain = Polygon(angles, vertices, parent)

    u1 = StandaloneVertexEigenfunction(vertex(domain, 1), fmpq(1//2), fmpq(3//2))
    u2 = StandaloneVertexEigenfunction(vertex(domain, 6), fmpq(0), fmpq(3//2))
    u3 = StandaloneVertexEigenfunction(vertex(domain, 7), fmpq(-1//2), fmpq(3//2))
    u4 = StandaloneVertexEigenfunction(vertex(domain, 12), fmpq(1), fmpq(3//2))


    u = CombinedEigenfunction(
        domain,
        [u1, u2, u3, u4],
        us_to_boundary = [
            setdiff(boundaries(domain), (12, 1)),
            setdiff(boundaries(domain), (5, 6)),
            setdiff(boundaries(domain), (6, 7)),
            setdiff(boundaries(domain), (11, 12))
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
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cos(θ), sin(θ)) for θ in (2parent(π)/n).*(0:n-1)]
    domain = Polygon(angles, vertices, parent)

    if !lightning
        if linked
            # FIXME: Linked VertexEigenfunctions cannot have separate
            # boundaries they are active on
            us = [
                LinkedEigenfunction(
                    [
                        StandaloneVertexEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                        for i in boundaries(domain)
                    ]
                )
            ]
            us_to_boundary = fill(BitSet([1]), length(us))
        else
            us = [
                StandaloneVertexEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                for i in boundaries(domain)
            ]
            us_to_boundary = [
                setdiff(boundaries(domain), (i, mod1(i - 1, length(boundaries(domain)))))
                for i in eachindex(us)
            ]
        end
    else
        if linked
            us = [
                LinkedEigenfunction(
                    [
                        StandaloneLightningEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                        for i in boundaries(domain)
                    ],
                    ones(n),
                ),
                StandaloneInteriorEigenfunction(domain, stride = n),
            ]

            us_to_boundary = fill(BitSet([1]), length(us))
        else
            us = [[
                StandaloneLightningEigenfunction(vertex(domain, i), i*(1 - θ) + θ*1//2, θ)
                for i in boundaries(domain)
            ]; StandaloneInteriorEigenfunction(domain)]

            us_to_boundary = fill(boundaries(domain), length(us))
        end
    end

    u = CombinedEigenfunction(domain, us, us_to_boundary = us_to_boundary)

    return domain, u
end

function example_domain_triangle_in_triangle(parent = RealField(precision(BigFloat)))
    domain1 = Triangle(fmpq(1//3), fmpq(1//3), parent)
    domain2 = TransformedDomain(
        Triangle(fmpq(1//3), fmpq(1//3), parent),
        fmpq(0),
        parent(0.25),
        SVector(cospi(fmpq(1//6), parent), sinpi(fmpq(1//6), parent))/4
    )

    domain = IntersectedDomain(domain1, domain2)

    u1 = StandaloneLightningEigenfunction(domain1, 1)

    u2 = LinkedEigenfunction(
        [
            StandaloneLightningEigenfunction(domain1, i)
            for i in 2:3
        ]
    )

    u3 = StandaloneLightningEigenfunction(domain2, 1, l = parent(0.15), outside = true)

    u4 = LinkedEigenfunction(
        [
            StandaloneLightningEigenfunction(domain2, i, l = parent(0.15), outside = true)
            for i in 2:3
        ]
    )

    u5 = StandaloneInteriorEigenfunction(domain)


    us = [u1, u2, u3, u4, u5]
    #us = vcat([u1], u2.us, [u3], u4.us, [u5])
    orders = [1, 1, 1, 1, 1]
    us_to_boundary = fill(1:6, length(us))

    u = CombinedEigenfunction(
        domain,
        us,
        orders,
        us_to_boundary = us_to_boundary,
    )

    return domain, u
end

function example_domain_ngon_in_ngon(
    n1,
    n2,
    parent = RealField(precision(BigFloat));
    linked = true,
)
    θ1 = fmpq((n1 - 2)//n1)
    angles1 = fill(θ1, n1)
    vertices1 = [(cos(θ), sin(θ)) for θ in (2parent(π)/n1).*(0:n1-1)]
    domain1 = Polygon(angles1, vertices1, parent)

    θ2 = fmpq((n2 - 2)//n2)
    angles2 = fill(θ2, n2)
    vertices2 = [(cos(θ), sin(θ)) for θ in (2parent(π)/n2).*(0:n2-1)]
    domain2 = TransformedDomain(
        Polygon(angles2, vertices2, parent),
        fmpq(0),
        parent(0.5),
        SVector(parent(0), parent(0)),
    )

    domain = IntersectedDomain(domain1, domain2)

    if linked
        us = [
            LinkedEigenfunction(
                [StandaloneLightningEigenfunction(domain1, i) for i in boundaries(domain1)]
            ),
            LinkedEigenfunction(
                [
                    StandaloneLightningEigenfunction(domain2, i, outside = true, l = parent(0.4))
                    for i in boundaries(domain2)
                ]
            ),
            StandaloneInteriorEigenfunction(domain, stride = 4)
        ]

        u = CombinedEigenfunction(
            domain,
            us,
            [3, 3, 2],
            us_to_boundary = fill(BitSet([1, 5]), length(us)),
        )
    else
        us = vcat(
            [StandaloneLightningEigenfunction(domain1, i) for i in boundaries(domain1)],
            [
                StandaloneLightningEigenfunction(domain2, i, outside = true, l = parent(0.5))
                for i in boundaries(domain2)
            ],
            [StandaloneInteriorEigenfunction(domain)]
        )

        u = CombinedEigenfunction(
            domain,
            us,
            [3, 3, 3, 3, 3, 3, 3, 3, 2],
        )
    end

    return domain, u
end

"""
    example_domain_goal(parent = RealField(precision(BigFloat)))

Return the domain given by a hexagon with 6 triangles cut out as well
as a corresponding eigenfunction.

Computing the first few eigenvalues in Matlab using a finite element
method gives
```
47.1851
49.2548
49.2668
54.8239
54.8292
59.6856
61.8997
78.3523
87.4298
87.4398
```
though it's likely not accurate to more than 1 or 2 digits.

"""
function example_domain_goal(parent = RealField(precision(BigFloat)))
    # The main domain is a hexagon
    n = 6
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cospi(θ, parent), sinpi(θ, parent)) for θ in fmpq(2//n).*(0:n-1)]
    exterior = Polygon(angles, vertices, parent)

    interiors = [
        TransformedDomain(
            Triangle(fmpq(1//3), fmpq(1//3), parent),
            fmpq(i//3),
            parent(0.25),
            SVector(cospi(fmpq(i//3 + 1//6), parent), sinpi(fmpq(i//3 + 1//6), parent))/4,
        )
        for i in 0:5
    ]

    domain = IntersectedDomain(exterior, interiors)

    # Expansions from the vertices of the hexagon
    u1 = LinkedEigenfunction(
        [
            StandaloneLightningEigenfunction(
                vertex(exterior, i),
                fmpq(mod(Rational(1 - θ//2 + (i - 1)*(1 - θ)), 2)),
                θ
            )
            for i in boundaries(exterior)
        ]
    )

    # Expansions from the outer tip of the triangles
    u2 = LinkedEigenfunction(
        [
            StandaloneLightningEigenfunction(d, 1, l = parent(0.15), outside = true)
            for d in interiors
        ]
    )

    # Expansions from the inner tips of the triangles
    u3 = LinkedEigenfunction(
        [
            StandaloneLightningEigenfunction(d, i, l = parent(0.15), outside = true)
            for d in interiors, i in 2:3
        ][:]
    )

    # Expansion from the center
    u4 = StandaloneInteriorEigenfunction(domain, stride = 6)

    us = [u1, u2, u3, u4]
    us_to_boundary = fill(BitSet([1, 7, 8, 9]), length(us))

    u = CombinedEigenfunction(
        domain,
        us,
        us_to_boundary = us_to_boundary,
    )

    return domain, u
end

"""
    example_domain_goal_v2(parent = RealField(precision(BigFloat)))

Return the domain given by a hexagon with 6 polygon cut out as well
as a corresponding eigenfunction.

Computing the first few eigenvalues in Matlab using a finite element
method gives
```
26.3259
56.8914
56.9041
56.9041
```

"""
function example_domain_goal_v2(
    parent = RealField(precision(BigFloat));
    even = true,
    a = 0.18,
    b = 0.1,
    c = 0.2,
)
    # The main domain is a hexagon
    n = 6
    θ = fmpq((n - 2)//n)
    angles = fill(θ, n)
    vertices = [(cospi(θ, parent), sinpi(θ, parent)) for θ in fmpq(2//n).*(0:n-1)]
    exterior = Polygon(angles, vertices, parent)

    points = [
        SVector(parent(1), parent(0)),
        SVector(parent(1), parent(0)),
        SVector(cospi(fmpq(1//3), parent), sinpi(fmpq(1//3), parent)),
        SVector(cospi(fmpq(1//3), parent), sinpi(fmpq(1//3), parent)),
    ]

    points[1] *= a + b
    points[2] *= a
    points[3] *= a
    points[4] *= a + b

    points[2] += 2*sinpi(fmpq(1//3), parent)*b.*SVector(cospi(fmpq(1//6), parent), sinpi(fmpq(1//6), parent))
    points[3] += 2*sinpi(fmpq(1//3), parent)*b.*SVector(cospi(fmpq(1//6), parent), sinpi(fmpq(1//6), parent))

    # Compute angles for interior polygons given the vertices
    #interior_angles = [
    #    let a = LinearAlgebra.norm(points[mod1(i - 1, length(points))] - points[mod1(i + 1, length(points))]),
    #    b = LinearAlgebra.norm(points[mod1(i + 1, length(points))] - points[i]),
    #    c = LinearAlgebra.norm(points[mod1(i - 1, length(points))] - points[i])
    #    acos((b^2 + c^2 - a^2)/(2b*c))
    #    end
    #    for i in eachindex(points)
    #]
    interior_angles = [1//3, 2//3, 2//3, 1//3]

    interiors = [
        TransformedDomain(
            Polygon(interior_angles, points, parent),
            fmpq(i//3),
            parent(1),
            c.*SVector(cospi(fmpq(i//3 + 1//6), parent), sinpi(fmpq(i//3 + 1//6), parent)),
        )
        for i in 0:5
    ]

    domain = IntersectedDomain(exterior, interiors)

    # Expansions from the vertices of the hexagon
    u1 = LinkedEigenfunction(
        [
            #StandaloneLightningEigenfunction(
            #    vertex(exterior, i),
            #    fmpq(mod(Rational(1 - θ//2 + (i - 1)*(1 - θ)), 2)),
            #    θ,
            #    even = even,
            #)
            StandaloneVertexEigenfunction(
                vertex(exterior, i),
                i*(1 - θ) + θ*1//2,
                θ,
                2,
            )
            for i in boundaries(exterior)
        ]
    )

    # Expansions from the outer part of the interior polygons
    u2 = LinkedEigenfunction(
        [
            #StandaloneLightningEigenfunction(
            #    d,
            #    i,
            #    outside = true,
            #    l = parent(0.01),
            #)
            StandaloneLightningEigenfunction(
                vertex(d, i),
                ifelse(i == 2, fmpq(4//3), fmpq(5//3)) + fmpq(d.rotation),
                2 - d.original.angles[i],
                l = parent(0.08),
                even = even,
            )
            for d in interiors, i in [2, 3]
        ][:]
    )

    # Expansions from the inner part of the interior polygons
    u3 = LinkedEigenfunction(
        [
            #StandaloneLightningEigenfunction(
            #    d,
            #    i,
            #    outside = true,
            #    l = parent(0.02),
            #)
            StandaloneLightningEigenfunction(
                vertex(d, i),
                ifelse(i == 1, fmpq(2//3), fmpq(0)) + fmpq(d.rotation),
                2 - d.original.angles[i],
                l = parent(0.1),
                even = even,
            )
            for d in interiors, i in [1, 4]
        ][:]
    )

    # Expansion from the center
    u4 = StandaloneInteriorEigenfunction(SVector(parent(0), parent(0)), stride = 6; even)

    # Interior expansions closer to outer boundaries
    u5 = LinkedEigenfunction(
        [
            StandaloneInteriorEigenfunction(
                parent(3//4).*SVector(cospi(fmpq(i//3), parent), sinpi(fmpq(i//3), parent)),
                fmpq(i//3),
                stride = 6,
                even = true
            )
            for i in 0:5
        ]
    )

    us = [u1, u2, u3, u5]
    #us = [u1.us..., u2.us..., u3.us..., u4, u5.us...]
    us_to_boundary = fill(BitSet([1, 7, 8, 10]), length(us))

    u = CombinedEigenfunction(
        domain,
        us,
        ifelse(even, [2, 2, 2, 1], [3, 3, 3, 1]),
        us_to_boundary = us_to_boundary,
    )

    return domain, u
end
