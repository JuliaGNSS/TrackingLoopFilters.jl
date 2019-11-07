module TrackingLoopFilters

    using DocStringExtensions
    using LinearAlgebra
    using Unitful: Hz

    export
        propagate,
        get_filtered_output,
        filter_loop,
        FirstOrderLF,
        SecondOrderBilinearLF,
        SecondOrderBoxcarLF,
        ThirdOrderBilinearLF,
        ThirdOrderBoxcarLF,
        AbstractLoopFilter

    abstract type AbstractLoopFilter end

    include("FirstOrderLF.jl")
    include("SecondOrderLF.jl")
    include("ThirdOrderLF.jl")

    """
    $(SIGNATURES)

    Propagates the state of loop filter and returns the calculated output of the loop filter
    as well as the propagated loop filter state.
    """
    function filter_loop(state::T, δθ, Δt, bandwidth) where T <: AbstractLoopFilter
        next_state = propagate(state, δθ, Δt, bandwidth)
        output = get_filtered_output(next_state, δθ, Δt, bandwidth)
        output, next_state
    end

end
