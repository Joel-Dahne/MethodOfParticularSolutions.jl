"""
    orientation(domain::AbstractPlanarDomain, i::Integer; reversed = false)

Return the orientation of the vertex `i` of the domain.

That is, the angle when going in positive direction from the x-axis
until hitting the outside boundary of the domain at the vertex.

If `reversed = false` then instead compute the angle going in the
negative direction.

Depending on the domain this returns either a rational number, in
which case it should be multiplied by `π` to give the actual angle, or
the angle directly.
"""
orientation(domain::AbstractPlanarDomain, i::Integer; reversed = false)
