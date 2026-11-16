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

natural_frequency(state::AbstractThirdOrderLF, bandwidth) = bandwidth / 0.7845

# `(ω₀, ω₀_assist)` of the third order loop and the assisting second order loop.
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

A third order bilinear loop assisted by a second order loop, e.g. an FLL-assisted PLL
(Kaplan & Hegarty, Figure 5.21(b)). It takes a `Tuple` `(δθ_high, δθ_low)` of a phase
error, as [`ThirdOrderBilinearLF`](@ref), and a frequency error in the units of the
output, e.g. cycles and Hz.

The bandwidth is a `Tuple` `(bandwidth_high, bandwidth_low)` of the loops' noise
bandwidths: `ω₀ = bandwidth_high / 0.7845` and `ω₀_assist = bandwidth_low / 0.53`
(Kaplan & Hegarty, Table 5.6). `bandwidth_low = 0` turns the assistance off.

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderAssistedBilinearLF()
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

# From the states, `x1` being the frequency at the start of the first integration: a
# filter that has not been propagated yet has no integration time.
_unpropagated_integration_time(x1, x2) = _in_seconds(zero(x1) / oneunit(x2))
ThirdOrderBilinearLF(x1, x2) = ThirdOrderBilinearLF(x1, x2, _unpropagated_integration_time(x1, x2))
ThirdOrderAssistedBilinearLF(x1, x2) =
    ThirdOrderAssistedBilinearLF(x1, x2, _unpropagated_integration_time(x1, x2))
ThirdOrderBoxcarLF(x1, x2) = ThirdOrderBoxcarLF(x1, x2, _unpropagated_integration_time(x1, x2))

"""
    ThirdOrderBilinearLF()

Construct a third order bilinear loop filter with zero initial state.
"""
function ThirdOrderBilinearLF()
    ThirdOrderBilinearLF(0.0Hz, 0.0Hz^2)
end

"""
    ThirdOrderAssistedBilinearLF()

Construct a third order assisted bilinear loop filter with zero initial state.
"""
function ThirdOrderAssistedBilinearLF()
    ThirdOrderAssistedBilinearLF(0.0Hz, 0.0Hz^2)
end

"""
    ThirdOrderBoxcarLF()

Construct a third order boxcar loop filter with zero initial state.
"""
function ThirdOrderBoxcarLF()
    ThirdOrderBoxcarLF(0.0Hz, 0.0Hz^2)
end

# The frequency state is the frequency `Δt` (bilinear) or `Δt / 2` (boxcar) after the end
# of the latest integration of `Δt`.
_frequency_state_offset(state::Union{ThirdOrderBilinearLF,ThirdOrderAssistedBilinearLF}, Δt) = Δt
_frequency_state_offset(state::ThirdOrderBoxcarLF, Δt) = Δt / 2

# A change of the integration time, e.g. after bit synchronisation, does not shift the
# frequency estimate.
function at_integration_time(state::T, Δt) where T <: Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}
    shift = _frequency_state_offset(state, Δt) - _frequency_state_offset(state, state.Δt)
    T(state.x1 + shift * state.x2, state.x2, Δt)
end

# The frequency state integrates the boxcar output `x2 + Δt * ω₀³ * δθ` of the frequency
# rate integrator, which includes the current `δθ` (Kaplan & Hegarty, Figure 5.19(b)).
function propagate_state(state::ThirdOrderBoxcarLF, δθ, Δt, ω₀)
    x2_boxcar = state.x2 + Δt * ω₀^3 * δθ
    typeof(state)(state.x1 + Δt * x2_boxcar + 1.1 * Δt * ω₀^2 * δθ, x2_boxcar, state.Δt)
end

# The frequency state integrates the bilinear output `x2 + Δt / 2 * ω₀³ * δθ` of the
# frequency rate integrator (Kaplan & Hegarty, Figure 5.20(c)).
function propagate_state(state::ThirdOrderBilinearLF, δθ, Δt, ω₀)
    Δx2 = Δt * ω₀^3 * δθ
    x2_bilinear = state.x2 + Δx2 / 2
    typeof(state)(state.x1 + Δt * x2_bilinear + 1.1 * Δt * ω₀^2 * δθ, state.x2 + Δx2, state.Δt)
end

# As the bilinear filter, with the frequency path of the assisting loop (Kaplan & Hegarty,
# Figure 5.21(b)).
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

function filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, ω₀)
    state.x1 + Δt / 2 * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ
end

function filtered_output(state::ThirdOrderAssistedBilinearLF, δθ, Δt, natural_frequencies)
    ω₀, ω₀_assist = natural_frequencies
    state.x1 + Δt / 2 * state.x2 +
        (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ[1] +
        (sqrt(2) * ω₀_assist * Δt / 2 + ω₀_assist^2 * Δt^2 / 4) * δθ[2]
end

# Both boxcar integrators include the current `δθ` (Kaplan & Hegarty, Figure 5.19(b)).
function filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, ω₀)
    state.x1 + Δt * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt + ω₀^3 * Δt^2) * δθ
end

predict_frequency(state::Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}, τ) =
    _in_frequency_units(state.x1 + (τ - _frequency_state_offset(state, state.Δt)) * state.x2)

get_frequency_rate(state::Union{AbstractThirdOrderLF,AbstractThirdOrderAssistedLF}) =
    _in_frequency_rate_units(state.x2)

is_assisted(::Type{<:AbstractThirdOrderAssistedLF}) = true
