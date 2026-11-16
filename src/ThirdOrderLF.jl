"""
$(TYPEDEF)

Abstract base type for third order loop filters.
"""
abstract type AbstractThirdOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

Abstract base type for assisted third order loop filters.
"""
abstract type AbstractThirdOrderAssistedLF <: AbstractLoopFilter end

"""
$(SIGNATURES)

Natural frequency `ω₀ = bandwidth / 0.7845` of a third order loop filter with the
noise bandwidth `bandwidth` (Kaplan & Hegarty, Table 5.6).
"""
natural_frequency(state::AbstractThirdOrderLF, bandwidth) = bandwidth / 0.7845

"""
$(SIGNATURES)

Natural frequencies `(ω₀, ω₀_assist)` of the third order loop and the assisting
second order loop of an assisted third order filter, each from its own noise
bandwidth in `(bandwidth_high, bandwidth_low)` (Kaplan & Hegarty, Table 5.6).
"""
natural_frequency(state::AbstractThirdOrderAssistedLF, bandwidth::Tuple{Any,Any}) =
    (bandwidth[1] / 0.7845, bandwidth[2] / 0.53)

"""
$(TYPEDEF)

Third order bilinear loop filter.

Uses bilinear transformation for improved frequency response.
The natural frequency is `ω₀ = bandwidth / 0.7845` (Kaplan & Hegarty, Table 5.6).

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderBilinearLF()
lf = ThirdOrderBilinearLF(; frequency = 1000.0Hz, frequency_rate = 1.0Hz / s)
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct ThirdOrderBilinearLF{T1,T2,T3} <: AbstractThirdOrderLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
    "Integration time of the latest propagation, which the frequency state refers to"
    Δt::T3
end

"""
$(TYPEDEF)

Third order bilinear loop filter with second order assistance.

Combines a third order bilinear loop with a second order assisted loop
for improved tracking performance. Accepts a `Tuple` `(δθ_high, δθ_low)` of the
high and low order discriminator outputs: `δθ_high` is a phase error, as for
[`ThirdOrderBilinearLF`](@ref), and `δθ_low` a frequency error in the units of
the filter output, e.g. a PLL discriminator output in cycles and an FLL
discriminator output in Hz (Kaplan & Hegarty, Figure 5.21(b)).

The bandwidth is a `Tuple` `(bandwidth_high, bandwidth_low)` that sets each loop
from its own noise bandwidth: `ω₀ = bandwidth_high / 0.7845` and
`ω₀_assist = bandwidth_low / 0.53`, from the third and second order loops'
noise bandwidths `0.7845 · ω₀` and `0.53 · ω₀` (Kaplan & Hegarty, Table 5.6).
In a frequency-assisted phase-locked loop this sets the PLL and the FLL
bandwidths independently; `bandwidth_low = 0` turns the assistance off.

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderAssistedBilinearLF()
lf = ThirdOrderAssistedBilinearLF(; frequency = 1000.0Hz, frequency_rate = 1.0Hz / s)
output, next_lf = filter_loop(lf, (δθ_high, δθ_low), Δt, (bandwidth_high, bandwidth_low))
```
"""
struct ThirdOrderAssistedBilinearLF{T1,T2,T3} <: AbstractThirdOrderAssistedLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
    "Integration time of the latest propagation, which the frequency state refers to"
    Δt::T3
end

"""
$(TYPEDEF)

Third order boxcar loop filter.

Uses boxcar (rectangular) integration for simpler implementation.
The natural frequency is `ω₀ = bandwidth / 0.7845` (Kaplan & Hegarty, Table 5.6).

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderBoxcarLF()
lf = ThirdOrderBoxcarLF(; frequency = 1000.0Hz, frequency_rate = 1.0Hz / s)
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct ThirdOrderBoxcarLF{T1,T2,T3} <: AbstractThirdOrderLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
    "Integration time of the latest propagation, which the frequency state refers to"
    Δt::T3
end

# A filter that has not been propagated yet has no integration time, so its frequency
# state is the frequency at the start of its first integration.
_unpropagated_integration_time(x1, x2) = _in_seconds(zero(x1) / oneunit(x2))

"""
    ThirdOrderBilinearLF(x1, x2)

Construct a third order bilinear loop filter with the frequency state `x1`, the
frequency at the start of its first integration, and the frequency rate state `x2`.
"""
ThirdOrderBilinearLF(x1, x2) = ThirdOrderBilinearLF(x1, x2, _unpropagated_integration_time(x1, x2))

"""
    ThirdOrderAssistedBilinearLF(x1, x2)

Construct a third order assisted bilinear loop filter with the frequency state `x1`,
the frequency at the start of its first integration, and the frequency rate state `x2`.
"""
ThirdOrderAssistedBilinearLF(x1, x2) =
    ThirdOrderAssistedBilinearLF(x1, x2, _unpropagated_integration_time(x1, x2))

"""
    ThirdOrderBoxcarLF(x1, x2)

Construct a third order boxcar loop filter with the frequency state `x1`, the
frequency at the start of its first integration, and the frequency rate state `x2`.
"""
ThirdOrderBoxcarLF(x1, x2) = ThirdOrderBoxcarLF(x1, x2, _unpropagated_integration_time(x1, x2))

"""
    ThirdOrderBilinearLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s)

Construct a third order bilinear loop filter at the initial `frequency` and
`frequency_rate`, `frequency` being the frequency at the start of the first
integration (see [`predict_frequency`](@ref)).
"""
ThirdOrderBilinearLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s) =
    ThirdOrderBilinearLF(_as_frequency_state(frequency), _as_frequency_rate_state(frequency_rate))

"""
    ThirdOrderAssistedBilinearLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s)

Construct a third order assisted bilinear loop filter at the initial `frequency` and
`frequency_rate`, `frequency` being the frequency at the start of the first
integration (see [`predict_frequency`](@ref)).
"""
ThirdOrderAssistedBilinearLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s) =
    ThirdOrderAssistedBilinearLF(_as_frequency_state(frequency), _as_frequency_rate_state(frequency_rate))

"""
    ThirdOrderBoxcarLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s)

Construct a third order boxcar loop filter at the initial `frequency` and
`frequency_rate`, `frequency` being the frequency at the start of the first
integration (see [`predict_frequency`](@ref)).
"""
ThirdOrderBoxcarLF(; frequency = 0.0Hz, frequency_rate = 0.0Hz / s) =
    ThirdOrderBoxcarLF(_as_frequency_state(frequency), _as_frequency_rate_state(frequency_rate))

# The frequency state of the bilinear filters is the frequency `Δt` after the end of
# the latest integration, that of the boxcar filter the frequency `Δt / 2` after it,
# `Δt` being the integration time of that integration.
_frequency_state_offset(state::Union{ThirdOrderBilinearLF,ThirdOrderAssistedBilinearLF}, Δt) = Δt
_frequency_state_offset(state::ThirdOrderBoxcarLF, Δt) = Δt / 2

"""
$(SIGNATURES)

The third order loop filter `state` with its frequency state moved to the moment it
refers to under the integration time `Δt`, which is `Δt` after the end of the latest
integration for the bilinear filters and `Δt / 2` for the boxcar filter. The state
stores the integration time it refers to, so a change of the integration time, e.g.
after bit synchronisation, does not shift its frequency estimate.
"""
function at_integration_time(state::T, Δt) where T <: Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}
    shift = _frequency_state_offset(state, Δt) - _frequency_state_offset(state, state.Δt)
    T(state.x1 + shift * state.x2, state.x2, Δt)
end

"""
$(SIGNATURES)

Propagate the third order boxcar loop filter state.

Updates both frequency and frequency rate state estimates. The frequency state
integrates the boxcar output `x2 + Δt * ω₀³ * δθ` of the frequency rate
integrator, which includes the current `δθ` (Kaplan & Hegarty, Figure 5.19(b)).

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
New loop filter state with updated estimates.
"""
function propagate_state(state::ThirdOrderBoxcarLF, δθ, Δt, ω₀)
    x2_boxcar = state.x2 + Δt * ω₀^3 * δθ
    typeof(state)(state.x1 + Δt * x2_boxcar + 1.1 * Δt * ω₀^2 * δθ, x2_boxcar, state.Δt)
end

"""
$(SIGNATURES)

Propagate the third order bilinear loop filter state.

Updates both frequency and frequency rate state estimates. The frequency state
integrates the bilinear output `x2 + Δt / 2 * ω₀³ * δθ` of the frequency rate
integrator (Kaplan & Hegarty, Figure 5.20(c)).

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
New loop filter state with updated estimates.
"""
function propagate_state(state::ThirdOrderBilinearLF, δθ, Δt, ω₀)
    Δx2 = Δt * ω₀^3 * δθ
    x2_bilinear = state.x2 + Δx2 / 2
    typeof(state)(state.x1 + Δt * x2_bilinear + 1.1 * Δt * ω₀^2 * δθ, state.x2 + Δx2, state.Δt)
end

"""
$(SIGNATURES)

Propagate the assisted third order loop filter state.

Updates both frequency and frequency rate state estimates using
dual discriminator inputs. The frequency state integrates the bilinear output
of the frequency rate integrator (Kaplan & Hegarty, Figure 5.21(b)).

# Arguments
- `state`: Current loop filter state
- `δθ`: `Tuple` `(δθ_high, δθ_low)` of the phase and frequency discriminator outputs
  (see [`ThirdOrderAssistedBilinearLF`](@ref))
- `Δt`: Integration time
- `natural_frequencies`: `Tuple` `(ω₀, ω₀_assist)` of the natural frequencies of the
  third order loop and the assisting second order loop

# Returns
New loop filter state with updated estimates.
"""
function propagate_state(state::T, δθ, Δt, natural_frequencies) where T <: AbstractThirdOrderAssistedLF
    ω₀, ω₀_assist = natural_frequencies
    Δx2 = Δt * ω₀^3 * δθ[1] + Δt * ω₀_assist^2 * δθ[2]
    x2_bilinear = state.x2 + Δx2 / 2
    T(
        state.x1 + Δt * x2_bilinear + 1.1 * Δt * ω₀^2 * δθ[1] + Δt * sqrt(2) * ω₀_assist * δθ[2],
        state.x2 + Δx2,
        state.Δt,
    )
end

"""
$(SIGNATURES)

Calculate the filtered output for the third order bilinear loop filter.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
Filtered frequency estimate.
"""
function filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, ω₀)
    state.x1 + Δt / 2 * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ
end

"""
$(SIGNATURES)

Calculate the filtered output for the assisted third order bilinear loop filter.

Combines outputs from the third order loop and second order assisted loop.

# Arguments
- `state`: Current loop filter state
- `δθ`: `Tuple` `(δθ_high, δθ_low)` of the phase and frequency discriminator outputs
  (see [`ThirdOrderAssistedBilinearLF`](@ref))
- `Δt`: Integration time
- `natural_frequencies`: `Tuple` `(ω₀, ω₀_assist)` of the natural frequencies of the
  third order loop and the assisting second order loop

# Returns
Filtered frequency estimate.
"""
function filtered_output(state::ThirdOrderAssistedBilinearLF, δθ, Δt, natural_frequencies)
    ω₀, ω₀_assist = natural_frequencies
    state.x1 + Δt / 2 * state.x2 +
        (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ[1] +
        (sqrt(2) * ω₀_assist * Δt / 2 + ω₀_assist^2 * Δt^2 / 4) * δθ[2]
end

"""
$(SIGNATURES)

Calculate the filtered output for the third order boxcar loop filter.

Both integrators include the current `δθ` (Kaplan & Hegarty, Figure 5.19(b)).
Calculate the output before [`propagate_state`](@ref) propagates the state with the
same `δθ`, as [`filter_loop`](@ref) does. From the propagated state, the current `δθ`
would be integrated twice.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
Filtered frequency estimate.
"""
function filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, ω₀)
    state.x1 + Δt * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt + ω₀^3 * Δt^2) * δθ
end

predict_frequency(state::Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}, τ) =
    _in_frequency_units(state.x1 + (τ - _frequency_state_offset(state, state.Δt)) * state.x2)

get_frequency_rate(state::Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}) =
    _in_frequency_rate_units(state.x2)

is_assisted(::Type{<:AbstractThirdOrderAssistedLF}) = true

_estimates_frequency_rate(::Type{<:Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}}) = true
