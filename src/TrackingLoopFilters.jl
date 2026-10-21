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

    Propagate the loop filter state.

    Calculate the output with [`get_filtered_output`](@ref) before propagating the state
    with the same `δθ`, or use [`filter_loop`](@ref), which does both.

    # Arguments
    - `state`: Current loop filter state
    - `δθ`: Phase discriminator output, or for the assisted filter a `Tuple` `(δθ_high, δθ_low)`
    - `Δt`: Integration time
    - `bandwidth`: Loop bandwidth, or for the assisted filter a `Tuple`
      `(bandwidth_high, bandwidth_low)` (see [`ThirdOrderAssistedBilinearLF`](@ref))

    # Returns
    New loop filter state.
    """
    function propagate(state::AbstractLoopFilter, δθ, Δt, bandwidth)
        propagate_state(state, δθ, Δt, natural_frequency(state, bandwidth))
    end

    """
    $(SIGNATURES)

    Calculate the filtered output for the current loop filter state.

    Calculate it before [`propagate`](@ref) propagates the state with the same `δθ`, as
    [`filter_loop`](@ref) does. From the propagated state, the current `δθ` would be
    integrated twice.

    # Arguments
    - `state`: Current loop filter state
    - `δθ`: Phase discriminator output, or for the assisted filter a `Tuple` `(δθ_high, δθ_low)`
    - `Δt`: Integration time
    - `bandwidth`: Loop bandwidth, or for the assisted filter a `Tuple`
      `(bandwidth_high, bandwidth_low)` (see [`ThirdOrderAssistedBilinearLF`](@ref))

    # Returns
    Filtered frequency estimate.
    """
    function get_filtered_output(state::AbstractLoopFilter, δθ, Δt, bandwidth)
        filtered_output(state, δθ, Δt, natural_frequency(state, bandwidth))
    end

    """
    $(SIGNATURES)

    Propagate the loop filter state and return both the filtered output and next state.

    This is a convenience function that combines [`propagate`](@ref) and
    [`get_filtered_output`](@ref) into a single call. The output is calculated from
    the current state, before it is propagated with the same `δθ`.

    # Arguments
    - `state`: Current loop filter state
    - `δθ`: Phase discriminator output, or for the assisted filter a `Tuple` `(δθ_high, δθ_low)`
    - `Δt`: Integration time
    - `bandwidth`: Loop bandwidth, or for the assisted filter a `Tuple`
      `(bandwidth_high, bandwidth_low)` (see [`ThirdOrderAssistedBilinearLF`](@ref))

    # Returns
    - `output`: Filtered frequency estimate
    - `next_state`: Propagated loop filter state
    """
    function filter_loop(state::T, δθ, Δt, bandwidth) where T <: AbstractLoopFilter
        ω₀ = natural_frequency(state, bandwidth)
        output = filtered_output(state, δθ, Δt, ω₀)
        next_state = propagate_state(state, δθ, Δt, ω₀)
        output, next_state
    end

end
