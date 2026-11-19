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
- [`filter_loop`](@ref): The filtered output and the propagated state
- [`predict_frequency`](@ref): The filter's frequency estimate at a given time
- [`get_frequency_rate`](@ref): The filter's frequency rate estimate
- [`reinit`](@ref): A filter of the same type at a given frequency and frequency rate
- [`is_assisted`](@ref): Whether the filter takes a second, frequency discriminator input
"""
module TrackingLoopFilters

    using DocStringExtensions
    using LinearAlgebra
    using Unitful: Hz, s, uconvert

    export
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

    # Frequencies in Hz and frequency rates in Hz/s, whatever units the filter
    # equations promoted the states to.
    _in_frequency_units(f) = uconvert(Hz, f)
    _in_frequency_rate_units(r) = uconvert(Hz / s, r)
    _in_seconds(t) = uconvert(s, t)

    # The filter states hold frequencies in Hz and frequency rates in Hz², the units
    # the filter equations leave them in, so that a filter keeps its type when it is
    # propagated.
    _as_frequency_state(f) = uconvert(Hz, float(f))
    _as_frequency_rate_state(r) = uconvert(Hz^2, float(r))

    """
        reinit(state_or_type; frequency = 0.0Hz, frequency_rate = 0.0Hz / s)

    A loop filter of the type of `state`, or of `type`, at the initial `frequency`
    and `frequency_rate`, `frequency` being the frequency at the start of the first
    integration (see [`predict_frequency`](@ref)). Of a concrete type, e.g. that of
    `state`, it keeps the element types and units. Filters below third order
    estimate no frequency rate and ignore `frequency_rate`, as
    [`get_frequency_rate`](@ref) returns zero for them, so generic code can pass a
    rate to any filter.

    To move the estimates of a filter, e.g. onto another reference frequency, pass
    its own estimates:

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

    # Whether filters of the type estimate a frequency rate.
    _estimates_frequency_rate(::Type{<:AbstractLoopFilter}) = false

    # The filter `lf` with its fields converted to those of the concrete type `F`.
    _with_field_types_of(::Type{F}, lf) where F = F(ntuple(i -> getfield(lf, i), Val(fieldcount(F)))...)

    """
        predict_frequency(state, τ)

    The loop filter's estimate of the frequency `τ` after the end of the latest
    integration it was propagated with (for a filter that has not been propagated
    yet: after the start of its first integration), extrapolated with its frequency
    rate estimate. `τ` is a time, e.g. `0.0s`, and may be negative, e.g. `-Δt / 2` for
    the centre of an integration of `Δt`. Requires states in units of frequency.

    The frequency states of different filters refer to different moments, which this
    accounts for: that of a bilinear third order filter is the frequency `Δt` after
    the end of the latest integration, that of the boxcar third order filter the
    frequency `Δt / 2` after it, `Δt` being the integration time of that integration,
    which the filter stores. Filters below third order estimate no frequency rate and
    return their frequency state for any `τ`.

    [`filter_loop`](@ref) calculates its output from the state before propagating it,
    so the output after the next integration commands the integration after it.
    Without a new discriminator output, and for integrations of the length `Δt` the
    filter was propagated with, that output is `predict_frequency(state, 3Δt / 2)`,
    the frequency at the centre of the integration it commands.
    """
    function predict_frequency end

    """
        get_frequency_rate(state)

    The loop filter's estimate of the frequency rate, in Hz/s. Zero for filters
    below third order, which estimate none. Requires states in units of frequency.
    """
    function get_frequency_rate end

    """
        is_assisted(state_or_type)

    Whether the loop filter is assisted by a second order frequency loop, and so
    takes a `Tuple` `(δθ_high, δθ_low)` of phase and frequency discriminator outputs
    (see [`ThirdOrderAssistedBilinearLF`](@ref)) instead of a single phase
    discriminator output.
    """
    is_assisted(::Type{<:AbstractLoopFilter}) = false
    is_assisted(state::AbstractLoopFilter) = is_assisted(typeof(state))

    """
    $(SIGNATURES)

    The loop filter `state` with its states moved to the moments they refer to under
    the integration time `Δt`. [`propagate`](@ref), [`get_filtered_output`](@ref) and
    [`filter_loop`](@ref) call it before `propagate_state` and `filtered_output`, which
    take a state at the integration time they are called with. Returns `state` for
    filters whose states refer to no moment that depends on the integration time.
    """
    at_integration_time(state::AbstractLoopFilter, Δt) = state

    include("FirstOrderLF.jl")
    include("SecondOrderLF.jl")
    include("ThirdOrderLF.jl")

    # The first and second order filters hold a constant frequency estimate. Adding their
    # zero frequency rate over `τ` raises a DimensionError for a `τ` that is not a time, as
    # for the third order filters.
    const ConstantFrequencyLF = Union{FirstOrderLF,AbstractSecondOrderLF}

    predict_frequency(state::ConstantFrequencyLF, τ) =
        _in_frequency_units(state.x + τ * get_frequency_rate(state))

    get_frequency_rate(state::ConstantFrequencyLF) = _in_frequency_rate_units(zero(state.x) / s)

    """
    $(SIGNATURES)

    Propagate the loop filter state. Not exported: use [`filter_loop`](@ref), which
    also calculates the output.

    Calculate the output with [`get_filtered_output`](@ref) before propagating the state
    with the same `δθ`.

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

    Calculate the filtered output for the current loop filter state. Not exported: use
    [`filter_loop`](@ref), which also propagates the state.

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

    The output is calculated from the current state, before it is propagated with the
    same `δθ`, so the output after the next integration commands the integration after
    it. A loop filter type defined outside the package defines `natural_frequency`,
    `propagate_state` and `filtered_output`, and `at_integration_time` if its states
    refer to a moment that depends on the integration time.

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

end
