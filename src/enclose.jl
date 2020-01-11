function maximize(u::AbstractEigenfunction,
                  λ::arb)
    @error "no rigorous implementation of maximize for $(typeof(u)), computing approximate maximum"
    boundary = boundary_points(u.domain, u, 1000)
    maximum(abs(u(b, λ)) for b in boundary)
end

function enclose_eigenvalue(domain::AbstractDomain,
                            u::AbstractEigenfunction,
                            λ::arb)
    m = maximize(u, λ)
    n = norm(u, λ)
    ϵ = sqrt(area(domain))*m/n

    enclosure = λ/ball(domain.parent(1), ϵ)
end
