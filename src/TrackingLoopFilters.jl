"""
    TrackingLoopFilters

A Julia package implementing loop filters for GNSS tracking algorithms.

Provides first, second, and third order loop filters with boxcar and bilinear
implementations. The filters take Unitful quantities: frequencies, e.g. in Hz, times,
e.g. in s, and bandwidths in units of frequency. Filtering without units is not
supported.

# Exported Types
- [`FirstOrderLF`](@ref): First order loop filter
- [`SecondOrderBilinearLF`](@ref): Second order bilinear loop filter
- [`SecondOrderBoxcarLF`](@ref): Second order boxcar loop filter
- [`ThirdOrderBilinearLF`](@ref): Third order bilinear loop filter
- [`ThirdOrderBoxcarLF`](@ref): Third order boxcar loop filter
- [`ThirdOrderAssistedBilinearLF`](@ref): Third order bilinear loop filter with second order assistance

# Exported Functions
- [`propagate`](@ref): Propagate the loop filter state
- [`get_filtered_output`](@ref): Get the filtered output for the current state
- [`filter_loop`](@ref): Combined propagate and get output in one call
- [`predict_frequency`](@ref): The filter's frequency estimate at a given time
- [`get_frequency_rate`](@ref): The filter's frequency rate estimate
- [`reinit`](@ref): A filter of the same type at a given frequency and frequency rate
- [`is_assisted`](@ref): Whether the filter takes a second, frequency discriminator input
"""
module TrackingLoopFilters

    using DocStringExtensions
    using Unitful: Hz, s, uconvert

    export
        propagate,
        get_filtered_output,
        filter_loop,
        predict_frequency,
        get_frequency_rate,
        is_assisted,
        reinit,
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

    # Estimates in Hz and Hz/s, whatever units the filter equations promoted the states to.
    _in_frequency_units(f) = uconvert(Hz, f)
    _in_frequency_rate_units(r) = uconvert(Hz / s, r)
    _in_seconds(t) = uconvert(s, t)

    # States in Hz and Hz², the units the filter equations leave them in, so that
    # propagating keeps the type of a filter.
    _as_frequency_state(f) = uconvert(Hz, float(f))
    _as_frequency_rate_state(r) = uconvert(Hz^2, float(r))

    """
        reinit(state_or_type; frequency = 0.0Hz, frequency_rate = 0.0Hz / s)

    A loop filter of the type of `state` or `type` at `frequency`, the frequency at the
    start of its first integration, and `frequency_rate`. A given filter or concrete type
    keeps its element types and units. Filters below third order ignore `frequency_rate`.

    To move the estimates of a filter, e.g. onto another reference frequency:

    ```julia
    reinit(lf; frequency = predict_frequency(lf, 0.0s) - Δf,
        frequency_rate = get_frequency_rate(lf) - Δfr)
    ```
    """
    function reinit end

    reinit(state::AbstractLoopFilter; kwargs...) = reinit(typeof(state); kwargs...)

    function reinit(::Type{F}; frequency = 0.0Hz, frequency_rate = 0.0Hz / s) where F <: AbstractLoopFilter
        G = Base.typename(F).wrapper
        lf = _estimates_frequency_rate(F) ? G(; frequency, frequency_rate) : G(; frequency)
        isconcretetype(F) ? _with_field_types_of(F, lf) : lf
    end

    _estimates_frequency_rate(::Type{<:AbstractLoopFilter}) = false

    _with_field_types_of(::Type{F}, lf) where F = F(ntuple(i -> getfield(lf, i), Val(fieldcount(F)))...)

    """
        predict_frequency(state, τ)

    The loop filter's frequency estimate `τ` after the end of its latest integration (before
    the first one: after its start), extrapolated with its frequency rate estimate. `τ` is a
    time and may be negative, e.g. `-Δt / 2` for the centre of an integration of `Δt`.

    Without a new discriminator output, the output of [`filter_loop`](@ref) after the next
    integration of `Δt` is `predict_frequency(state, 3Δt / 2)`, the frequency at the centre
    of the integration it commands.
    """
    function predict_frequency end

    """
        get_frequency_rate(state)

    The loop filter's frequency rate estimate in Hz/s; zero below third order.
    """
    function get_frequency_rate end

    """
        is_assisted(state_or_type)

    Whether the loop filter is assisted by a second order frequency loop and so takes a
    `Tuple` `(δθ_high, δθ_low)` of discriminator outputs (see
    [`ThirdOrderAssistedBilinearLF`](@ref)).
    """
    is_assisted(::Type{<:AbstractLoopFilter}) = false
    is_assisted(state::AbstractLoopFilter) = is_assisted(typeof(state))

    """
        at_integration_time(state, Δt)

    The loop filter `state` with its frequency state moved to the moment it refers to for
    the integration time `Δt`, before `propagate_state` and `filtered_output`. Returns
    `state` for filters whose states do not depend on the integration time.
    """
    at_integration_time(state::AbstractLoopFilter, Δt) = state

    include("FirstOrderLF.jl")
    include("SecondOrderLF.jl")
    include("ThirdOrderLF.jl")

    # No frequency rate state: the prediction holds the frequency over `τ`. Adding the zero
    # rate over `τ` raises a DimensionError for a `τ` that is not a time, as for the third
    # order filters.
    const NoFrequencyRateLF = Union{FirstOrderLF,AbstractSecondOrderLF}

    predict_frequency(state::NoFrequencyRateLF, τ) =
        _in_frequency_units(state.x + τ * get_frequency_rate(state))

    get_frequency_rate(state::NoFrequencyRateLF) = _in_frequency_rate_units(zero(state.x) / s)

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
        propagate_state(at_integration_time(state, Δt), δθ, Δt, natural_frequency(state, bandwidth))
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
        filtered_output(at_integration_time(state, Δt), δθ, Δt, natural_frequency(state, bandwidth))
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
        state = at_integration_time(state, Δt)
        ω₀ = natural_frequency(state, bandwidth)
        output = filtered_output(state, δθ, Δt, ω₀)
        next_state = propagate_state(state, δθ, Δt, ω₀)
        output, next_state
    end

    """
        natural_frequency(state, bandwidth)

    The natural frequency `ω₀` of the loop filter for the noise `bandwidth` (Kaplan &
    Hegarty, Table 5.6).
    """
    function natural_frequency end

    """
        propagate_state(state, δθ, Δt, ω₀)

    The loop filter state propagated with the discriminator output `δθ` of an integration
    of `Δt`, for the natural frequency `ω₀`.
    """
    function propagate_state end

    """
        filtered_output(state, δθ, Δt, ω₀)

    The loop filter output for the discriminator output `δθ` of an integration of `Δt`,
    calculated from the state before `propagate_state` propagates it with the same `δθ`,
    which would otherwise integrate `δθ` twice.
    """
    function filtered_output end

end
