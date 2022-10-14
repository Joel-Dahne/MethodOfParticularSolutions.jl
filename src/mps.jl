"""
    mps!(u, domain, a, b, N; num_boundary = 2N, num_interior = 2N, optim_prec = 20)

Return the `λ` minimizing `σ(λ)` on the interval `[a, b]` and set `u`
to the corresponding eigenfunction.

The number of coefficients used for `u` is determined by `N`. The
number of boundary points and interior points used in the calculations
are given by `num_boundary` and `num_interior` respectively.

The minimum is computer to `optim_prec` bits of precision.

If `show_progress` is true then show a progress bar for the
minimization. If `show_trace` is true then show the trace of the
minimization. If `store_trace` is true also return the result of the
search for the minimum, if `a == b` (in which case no search is
performed) return `nothing`.
"""
function mps!(
    u::AbstractEigenfunction,
    domain::AbstractDomain,
    a::Real,
    b::Real,
    N::Integer = 8;
    num_boundary::Integer = 2N,
    num_interior::Integer = 2N,
    optim_prec::Integer = 20,
    show_progress::Bool = false,
    show_trace::Bool = false,
    store_trace::Bool = false,
)
    # Some types of eigenfunctions (e.g.
    # StandaloneLightningEigenfunction) needs that the number of terms
    # used is determined from the beginning.
    set_eigenfunction!(u, zeros(N))

    a > b && throw(ArgumentError("must have a <= b, got a = $a, b = $b"))

    if a != b
        # Compute minimum of σ(λ)
        σ(λ) = sigma(λ, domain, u, N; num_boundary, num_interior)

        tol = 2.0^(-(optim_prec - 1))

        if show_progress
            start_prec = -log2((b - a) / ((b + a) / 2))

            optim_progress(state) = begin
                if state isa AbstractArray
                    state = state[end]
                end
                if state == :done
                    progress = "done"
                else
                    abs_error = state.metadata["x_upper"] - state.metadata["x_lower"]
                    rel_prec = -log2(abs_error / state.metadata["minimizer"])
                    progress = Float64(
                        max((rel_prec - start_prec) / (optim_prec - start_prec), 0),
                    )
                end
                @info "Computing minimum of σ(λ)" progress = progress
                false
            end
        else
            optim_progress = nothing
        end

        res = optimize(
            σ,
            a,
            b,
            rel_tol = tol,
            abs_tol = tol,
            callback = optim_progress;
            show_trace,
            store_trace,
        )

        if show_progress
            optim_progress(:done)
        end

        if !Optim.converged(res)
            @warn "Failed to compute minimum of σ(λ)"
        end

        λ = res.minimizer
    else
        λ = a
    end

    # Compute the eigenfunction corresponding to the computed minimum
    # of σ(λ).
    coefficients = sigma_coefficients(λ, domain, u, N; num_boundary, num_interior)
    set_eigenfunction!(u, coefficients)

    if store_trace
        if a != b
            return λ, res
        else
            return λ, nothing
        end
    else
        return λ
    end
end

"""
    iteratemps(domain, u, enclosure, Ns; kwargs...)

For each `N` in `Ns` compute the minimizing `λ` in `enclosure` and
corresponding eigenfunction with [`mps`](@ref). At each step also
compute an enclosure for `λ` and if this is better than the original
enclosure use that in the next step instead.

It takes a number of key-word arguments for handling progressively
updating parameters.
- `num_boundary_factor`, `num_interior_factor`: The number of points
  used on the boundary respectively in the interior is given by these
  factors multiplied by `N`.

There are a number of different ways for determining the tolerance to
use in the optimization of `σ(λ)` in the next step. It's
1. Constant: Use the tolerance every iteration. Set
`optim_prec_final` to the tolerance to use.
2. Linear: Specify a tolerance to use in the last iteration and then
linearly interpolate for the current `N`. Set `optim_prec_final` to
the tolerance and `optim_prec_linear` to true.
3. Step: Compute the current precision of the enclosure and use that
plus a fixed number as the tolerance for the next iteration. Set
`optim_prec_step` to true.
4. Adaptive: Estimate the improvement for the last step by computing
the precision of the current enclosure and the one in the iteration
before. Add the improvement to the precision of the current enclosure
plus a small constant. Set `optim_prec_adaptive` to true and
`optim_prec_adaptive_extra` to the extra precision to use.

The precision used in the computations is determined by the tolerance
for the optimization plus `extra_prec`.

If `recompute` is true then run `recompute!(u)` for each iteration.
This is useful for eigenfunction which precompute some of there
values, they can then be recomputed with the new precision. For most
eigenfunctions this doesn't have an effect.

If `rigorous_enclosure` is true then rigorously compute the enclosure
of the eigenvalue by rigorously bounding the maximum on the boundary
and the norm. If `rigorous_enclosure` is true then `rigorous_norm` can
be set to false to only compute the maximum rigorously.
"""
function iteratemps(
    domain::AbstractDomain,
    u::AbstractEigenfunction,
    enclosure::arb,
    Ns;
    num_boundary_factor::Integer = 2,
    num_interior_factor::Integer = 2,
    optim_prec_final = nothing,
    optim_prec_linear = false,
    optim_prec_step = nothing,
    optim_prec_adaptive = false,
    optim_prec_adaptive_extra = 10,
    extra_prec = 20,
    recompute = true,
    rigorous_enclosure = true,
    rigorous_norm = rigorous_enclosure,
    show_progress = false,
    show_trace = false,
    store_trace = false,
    extended_trace = false,
)
    # Make sure that at least one method for computing the precision
    # to use for the minimization is specified
    if isnothing(optim_prec_final) && isnothing(optim_prec_step) && !optim_prec_adaptive
        throw(
            ArgumentError(
                "at least one method for determining precision must be specified",
            ),
        )
    end

    # Setup function used for computing precision to which σ(λ) should
    # be minimized
    get_optim_prec(N, enclosure, λs, i) = begin
        optim_prec = typemax(Int)

        # If a final precision for the minimization is specified use
        # that, possibly with a linear interpolation.
        if !isnothing(optim_prec_final)
            if optim_prec_linear
                optim_prec = min(optim_prec, ceil(Int, optim_prec_final * N / Ns[end]))
            else
                optim_prec = min(optim_prec, optim_prec_final)
            end
        end

        # Set the precision to the relative accuracy of the current
        # enclosure plus a fixed precision
        if !isnothing(optim_prec_step)
            optim_prec =
                min(optim_prec, ArbTools.rel_accuracy_bits(enclosure) + optim_prec_step)
        end

        # Compute the precision increase between the last two
        # iterations and use the same plus a small constant
        if optim_prec_adaptive
            if i > 3 && isfinite(λs[i-2]) && isfinite(λs[i-1])
                optim_prec = min(
                    optim_prec,
                    max(
                        2ArbTools.rel_accuracy_bits(λs[i-1]) -
                        ArbTools.rel_accuracy_bits(λs[i-2]),
                        0,
                    ) + optim_prec_adaptive_extra,
                )
            else
                optim_prec = min(optim_prec, precision(domain.parent) - extra_prec)
            end
        end

        return optim_prec
    end


    # Copy the domain and eigenfunction to avoid changing the originals
    domain = deepcopy(domain)

    λs = similar(Ns, arb)
    us = [deepcopy(u) for _ in eachindex(Ns)]

    trace = MPSTrace()
    tracing = store_trace || show_trace || extended_trace

    if show_trace
        show(trace)
    end

    try
        for (i, N) in enumerate(Ns)
            ### Number of boundary and interior points to use ###
            num_boundary = num_boundary_factor * N
            num_interior = num_interior_factor * N

            ### Precision for optimization ###
            optim_prec = get_optim_prec(N, enclosure, λs, i)

            ### Precision to use in the computations ###
            new_prec = max(53, optim_prec + extra_prec) # Always use at least 53 bits

            ### Update the domain and eigenfunction with the new precision ###
            parent = RealField(new_prec)
            domain = typeof(domain)(domain; parent)
            set_domain!(us[i], domain)
            if recompute
                recompute!(us[i])
            end

            ### Run the computations ###
            @timeit_debug "$N" begin
                @timeit_debug "mps" begin
                    λ = setprecision(BigFloat, new_prec) do
                        mps!(
                            us[i],
                            domain,
                            getinterval(BigFloat, enclosure)...,
                            N,
                            store_trace = extended_trace;
                            num_boundary,
                            num_interior,
                            optim_prec,
                            show_progress,
                        )
                    end

                    if extended_trace
                        λ, optim_res = λ
                    end

                    λ = parent(λ)
                end

                @timeit_debug "enclosure" begin
                    if rigorous_enclosure
                        new_enclosure = enclose_eigenvalue(
                            domain,
                            us[i],
                            λ,
                            store_trace = tracing;
                            rigorous_norm,
                            extended_trace,
                            show_progress,
                        )
                    else
                        new_enclosure = enclose_eigenvalue_approx(
                            domain,
                            us[i],
                            λ,
                            max_numpoints = 4num_boundary_factor * N,
                            norm_numpoints = 4num_interior_factor * N,
                            store_trace = tracing;
                            extended_trace,
                        )
                    end
                end
            end

            if tracing
                metadata = Dict()

                if extended_trace
                    new_enclosure, n, m, maximize_trace = new_enclosure
                    metadata["optim_res"] = optim_res
                    metadata["maximize_trace"] = maximize_trace
                else
                    new_enclosure, n, m = new_enclosure
                end

                state = MPSState(N, new_prec, optim_prec, n, m, λ, new_enclosure, metadata)

                update!(trace, state, store_trace, show_trace)
            end

            λs[i] = new_enclosure
            if isfinite(new_enclosure)
                enclosure = setintersection(enclosure, new_enclosure)
            end
        end
    catch e
        if !(e isa InterruptException)
            rethrow(e)
        end
    end

    if store_trace
        return λs, us, trace
    else
        return λs, us
    end
end
