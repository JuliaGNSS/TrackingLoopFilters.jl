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

# Natural frequencies `(ω₀, ω₀_assist)` of the third order loop and the assisting
# second order loop of an assisted third order filter. A single noise bandwidth
# ties the assisting loop to the third order one (`ω₀ / 4`); a pair
# `(bandwidth_high, bandwidth_low)` sets each from its own noise bandwidth.
assisted_natural_frequencies(bandwidth) = (ω₀ = bandwidth * 1.2; (ω₀, ω₀ / 4))
assisted_natural_frequencies(bandwidth::Tuple{Any,Any}) = (bandwidth[1] * 1.2, bandwidth[2] / 0.53)

"""
$(TYPEDEF)

Third order bilinear loop filter.

Uses bilinear transformation for improved frequency response.
The natural frequency scaling factor is 1.2.

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderBilinearLF()
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct ThirdOrderBilinearLF{T1,T2} <: AbstractThirdOrderLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
end

"""
$(TYPEDEF)

Third order bilinear loop filter with second order assistance.

Combines a third order bilinear loop with a second order assisted loop
for improved tracking performance. Accepts a two-element discriminator
input vector `[δθ_high, δθ_low]`.

The bandwidth is either a single noise bandwidth, which sets the third order
loop and ties the assisting loop to it (`ω₀_assist = ω₀ / 4`), or a pair
`(bandwidth_high, bandwidth_low)` that sets each loop from its own noise
bandwidth: `ω₀ = 1.2 · bandwidth_high` and `ω₀_assist = bandwidth_low / 0.53`,
the second order loop's noise bandwidth `0.53 · ω₀` (Kaplan & Hegarty, Table
5.6). In a frequency-assisted phase-locked loop this sets the PLL and the FLL
bandwidths independently; `bandwidth_low = 0` turns the assistance off.

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderAssistedBilinearLF()
output, next_lf = filter_loop(lf, [δθ_high, δθ_low], Δt, bandwidth)
output, next_lf = filter_loop(lf, [δθ_high, δθ_low], Δt, (bandwidth_high, bandwidth_low))
```
"""
struct ThirdOrderAssistedBilinearLF{T1,T2} <: AbstractThirdOrderAssistedLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
end

"""
$(TYPEDEF)

Third order boxcar loop filter.

Uses boxcar (rectangular) integration for simpler implementation.
The natural frequency scaling factor is 1.2.

$(TYPEDFIELDS)

# Example
```julia
lf = ThirdOrderBoxcarLF()
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct ThirdOrderBoxcarLF{T1,T2} <: AbstractThirdOrderLF
    "Frequency state estimate"
    x1::T1
    "Frequency rate state estimate"
    x2::T2
end

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

"""
$(SIGNATURES)

Propagate the third order loop filter state.

Updates both frequency and frequency rate state estimates.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
New loop filter state with updated estimates.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractThirdOrderLF
    ω₀ = bandwidth * 1.2
    T(state.x1 + Δt * state.x2 + 1.1 * Δt * ω₀^2 * δθ, state.x2 + Δt * ω₀^3 * δθ)
end

"""
$(SIGNATURES)

Propagate the assisted third order loop filter state.

Updates both frequency and frequency rate state estimates using
dual discriminator inputs.

# Arguments
- `state`: Current loop filter state
- `δθ`: Two-element vector `[δθ_high, δθ_low]` with high and low order discriminator outputs
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth, or a pair `(bandwidth_high, bandwidth_low)` for the
  third order loop and the assisting second order loop (see
  [`ThirdOrderAssistedBilinearLF`](@ref))

# Returns
New loop filter state with updated estimates.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractThirdOrderAssistedLF
    ω₀, ω₀_assist = assisted_natural_frequencies(bandwidth)
    T(state.x1 + Δt * state.x2 + 1.1 * Δt * ω₀^2 * δθ[1] + Δt * sqrt(2) * ω₀_assist * δθ[2], state.x2 + Δt * ω₀^3 * δθ[1] + Δt * ω₀_assist^2 * δθ[2])
end

"""
$(SIGNATURES)

Calculate the filtered output for the third order bilinear loop filter.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀= bandwidth * 1.2
    state.x1 + Δt / 2 * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ
end

"""
$(SIGNATURES)

Calculate the filtered output for the assisted third order bilinear loop filter.

Combines outputs from the third order loop and second order assisted loop.

# Arguments
- `state`: Current loop filter state
- `δθ`: Two-element vector `[δθ_high, δθ_low]` with high and low order discriminator outputs
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth, or a pair `(bandwidth_high, bandwidth_low)` for the
  third order loop and the assisting second order loop (see
  [`ThirdOrderAssistedBilinearLF`](@ref))

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::ThirdOrderAssistedBilinearLF, δθ, Δt, bandwidth)
    ω₀, ω₀_assist = assisted_natural_frequencies(bandwidth)
    state.x1 + Δt / 2 * state.x2 +
        (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ[1] +
        (sqrt(2) * ω₀_assist * Δt / 2 + ω₀_assist^2 * Δt^2 / 4) * δθ[2]
end

"""
$(SIGNATURES)

Calculate the filtered output for the third order boxcar loop filter.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= bandwidth * 1.2
    state.x1 + 2.4 * ω₀ * δθ
end
