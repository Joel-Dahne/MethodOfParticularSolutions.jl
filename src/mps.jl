function mps(domain::AbstractDomain,
             eigenfunction::AbstractEigenfunction,
             enclosure::arb,
             N::Integer = 8;
             num_boundary = 2N,
             num_interior = 2N,
             optim_prec::Int = prec(domain.parent),
             norm_rigorous = true,
             store_trace = false,
             show_trace = false,
             show_progress = false,
             extended_trace = false,
             optim_store_trace = false,
             optim_show_trace = false,
             optim_extended_trace = false,
             enclose_store_trace = false,
             enclose_show_trace = false,
             enclose_extended_trace = false)
    # Compute minimum of σ(λ)
    σ = λ -> sigma(λ, domain, eigenfunction, N,
                   num_boundary = num_boundary,
                   num_interior = num_interior)

    tol = BigFloat(domain.parent(2)^(-(optim_prec - 1)))

    if show_progress
        start_prec = ArbTools.rel_accuracy_bits(enclosure)
        optim_progress(state) = begin
            if state == :done
                progress = "done"
            else
                abs_error = state.metadata["x_upper"] - state.metadata["x_lower"]
                rel_prec = -log2(abs_error/state.metadata["minimizer"])
                progress = max((rel_prec - start_prec)/(optim_prec - start_prec), 0)
            end
            @info "Computing minimum of σ(λ)" progress = progress
            false
        end
    else
        optim_progress = nothing
    end

    res = optimize(σ,
                   getinterval(BigFloat, enclosure)...,
                   rel_tol = tol,
                   abs_tol = tol,
                   store_trace = optim_store_trace,
                   show_trace = optim_show_trace,
                   extended_trace = optim_extended_trace,
                   callback = optim_progress)

    if show_progress
        optim_progress(:done)
    end

    if !Optim.converged(res)
        @warn "Failed to compute minimum of σ(λ)"
    end

    λ = res.minimizer

    # Compute the eigenfunction corresponding to the computed
    # minimum of σ(λ).
    coefficients = sigma_coefficients(λ, domain, eigenfunction, N,
                                      num_boundary = num_boundary,
                                      num_interior = num_interior)
    set_eigenfunction!(eigenfunction, coefficients)

    # Compute the enclosure of the eigenvalue
    λ = enclose_eigenvalue(domain,
                           eigenfunction,
                           domain.parent(λ),
                           norm_rigorous = norm_rigorous,
                           store_trace = enclose_store_trace,
                           show_trace = enclose_show_trace,
                           show_progress = show_progress,
                           extended_trace = enclose_extended_trace)

    return λ, eigenfunction
end

"""
    iteratemps(domain::AbstractDomain,
               eigenfunction::AbstractEigenfunction,
               enclosure::arb,
               Ns;
               kwargs...)
Call `mps` iteratively for the each value in `Ns`. Each time the
enclosure of the eigenvalue is refined.

It takes a number of key-word arguments for handling progressively
updating parameters.
- `num_boundary_factor`, `num_interior_factor`: The number of points
  used on the boundary respectively in the interior is given by these
  factors multiplied by `N`.

There are a number of different ways for determining the tolerance to
use in the optimization of σ(λ) in the next step. It's
1. Constant: Use the tolerance every iteration. Set
`optim_prec_final`.
2. Linear: Specify a tolerance to use in the last iteration and then
linearly interpolate for the current `N`. Set `optim_prec_final` and
`optim_prec_linear`.
3. Step: Compute the current precision of the enclosure and use that
plus a fixed number as the tolerance for the next iteration. Set
`optim_prec_step`.
4. Adaptive: Estimate the improvement for the last step by computing
the precision of the current enclosure and the one in the iteration
before. Add the improvement to the precision of the current enclosure
plus a small constant. Set `optim_prec_adaptive` and
`optim_prec_adaptive_extra`.

The precision used in the computations is determined by the tolerance
for the optimization plus `extra_prec`.
"""
function iteratemps(domain::AbstractDomain,
                    eigenfunction::AbstractEigenfunction,
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
                    norm_rigorous = true,
                    show_trace = false,
                    show_progress = false)
    # Make sure that at least one method for computing the precision
    # to use for the minimization is specified
    if isnothing(optim_prec_final) && isnothing(optim_prec_step) && !optim_prec_adaptive
        throw(ArgumentError("at least one method for determining precision must be specified"))
    end

    # Copy the domain and eigenfunction to avoid changing the originals
    domain = deepcopy(domain)
    eigenfunction = deepcopy(eigenfunction)

    λs = Vector{arb}(undef, length(Ns))

    if show_trace
        println("N    optim_prec    prec    λ    accuracy")
    end

    try
        for i in 1:length(Ns)
            N = Ns[i]

            ### Compute number of boundary and interior points to use ###
            num_boundary = num_boundary_factor*N
            num_interior = num_interior_factor*N

            ### Compute precision to which σ(λ) should be minimized ###
            optim_prec = typemax(Int)

            # If a final precision for the minimization is specified use
            # that, possibly with a linear interpolation.
            if !isnothing(optim_prec_final)
                if optim_prec_linear
                    optim_prec = min(optim_prec,
                                     ceil(Int, optim_prec_final*N/Ns[end]))
                else
                    optim_prec = min(optim_prec, optim_prec_final)
                end
            end

            # Set the precision to the relative accuracy of the current
            # enclosure plus a fixed precision
            if !isnothing(optim_prec_step)
                optim_prec = min(optim_prec,
                                 ArbTools.rel_accuracy_bits(enclosure) + optim_prec_step)
            end

            # Compute the precision increase between the last two
            # iterations and use the same plus a small constant
            if optim_prec_adaptive
                if i > 3 && isfinite(λs[i - 2]) && isfinite(λs[i - 1])
                    optim_prec = min(optim_prec,
                                     max(2ArbTools.rel_accuracy_bits(λs[i - 1])
                                         - ArbTools.rel_accuracy_bits(λs[i - 2]),
                                         0)
                                     + optim_prec_adaptive_extra)
                else
                    optim_prec = min(optim_prec,
                                     prec(domain.parent) - extra_prec)
                end
            end

            ### Compute precision to use in the computations ###
            # Always use at least 53 bits of precision
            new_prec = max(53, optim_prec + extra_prec)
            RR = RealField(new_prec)
            domain = typeof(domain)(domain, RR)
            eigenfunction.domain = typeof(eigenfunction.domain)(eigenfunction.domain, RR)

            ### Run the computations ###
            @timeit_debug "$N" begin
                λ, _ = setprecision(BigFloat, new_prec) do
                    mps(domain, eigenfunction, enclosure, N,
                        num_boundary = num_boundary,
                        num_interior = num_interior,
                        optim_prec = optim_prec,
                        norm_rigorous = norm_rigorous,
                        show_progress = show_progress)
                end
            end

            λs[i] = λ
            if isfinite(λ)
                enclosure = setintersection(enclosure, λ)
            end

            if show_trace
                println("$N    $optim_prec    $(prec(domain.parent))    $λ    $(ArbTools.rel_accuracy_bits(λ))")
            end
        end
    catch e
        if e isa InterruptException
            return λs
        else
            rethrow(e)
        end
    end

    return λs
end
