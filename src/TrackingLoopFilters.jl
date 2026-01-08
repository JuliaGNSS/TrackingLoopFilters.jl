"""
    TrackingLoopFilters

A Julia package implementing loop filters for GNSS tracking algorithms.

Provides first, second, and third order loop filters with boxcar and bilinear
implementations. All filters support Unitful quantities for type-safe calculations.

# Exported Types
- [`FirstOrderLF`](@ref): First order loop filter (stateless)
- [`SecondOrderBilinearLF`](@ref): Second order bilinear loop filter
- [`SecondOrderBoxcarLF`](@ref): Second order boxcar loop filter
- [`ThirdOrderBilinearLF`](@ref): Third order bilinear loop filter
- [`ThirdOrderBoxcarLF`](@ref): Third order boxcar loop filter
- [`ThirdOrderAssistedBilinearLF`](@ref): Third order bilinear loop filter with second order assistance

# Exported Functions
- [`propagate`](@ref): Propagate the loop filter state
- [`get_filtered_output`](@ref): Get the filtered output for the current state
- [`filter_loop`](@ref): Combined propagate and get output in one call
"""
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
        ThirdOrderAssistedBilinearLF,
        ThirdOrderBoxcarLF,
        AbstractLoopFilter

    """
    $(TYPEDEF)

    Abstract base type for all loop filters.
    """
    abstract type AbstractLoopFilter end

    include("FirstOrderLF.jl")
    include("SecondOrderLF.jl")
    include("ThirdOrderLF.jl")

    """
    $(SIGNATURES)

    Propagate the loop filter state and return both the filtered output and next state.

    This is a convenience function that combines [`propagate`](@ref) and
    [`get_filtered_output`](@ref) into a single call.

    # Arguments
    - `state`: Current loop filter state
    - `δθ`: Phase discriminator output (dimensionless for standard filters, vector for assisted)
    - `Δt`: Integration time
    - `bandwidth`: Loop bandwidth

    # Returns
    - `output`: Filtered frequency estimate
    - `next_state`: Propagated loop filter state
    """
    function filter_loop(state::T, δθ, Δt, bandwidth) where T <: AbstractLoopFilter
        next_state = propagate(state, δθ, Δt, bandwidth)
        output = get_filtered_output(next_state, δθ, Δt, bandwidth)
        output, next_state
    end

end
