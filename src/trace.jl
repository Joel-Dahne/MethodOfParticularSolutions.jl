struct MPSState
    N::Int
    precision::Int
    optim_prec::Int
    norm::arb
    maximum::arb
    minimizer::arb
    enclosure::arb
    metadata::Dict
end

struct MPSTrace
    states::Vector{MPSState}
end

MPSTrace() = MPSTrace(MPSState[])

function Base.show(io::IO, st::MPSState)
    @printf io "%4d    %4d    %9d    %8d    " st.N st.precision st.optim_prec ArbTools.rel_accuracy_bits(
        st.enclosure,
    )
    @printf io "%.5f    %26s    %s\n" Float64(st.norm) ArbTools.format_arb(
        st.maximum / st.norm,
        5,
    ) string(st.enclosure)
    return
end

function Base.show(io::IO, tr::MPSTrace)
    @printf io "   N    Prec     Opt prec    Enc prec    "
    @printf io "   Norm                  Maximum/Norm    Enclosure\n"
    @printf io "----    ----     --------    --------    "
    @printf io "-------    --------------------------    ---------\n"
    for st in tr.states
        show(io, st)
    end
    return
end

Base.push!(tr::MPSTrace, st::MPSState) = push!(tr.states, st)

function update!(
    tr::MPSTrace,
    st::MPSState,
    store_trace::Bool,
    show_trace::Bool,
    callback = nothing,
)
    if store_trace
        push!(tr, st)
    end
    if show_trace
        show(st)
    end
    if !isnothing(callback)
        callback(st)
    end
    return
end
