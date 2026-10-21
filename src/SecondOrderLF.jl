"""
$(TYPEDEF)

Abstract base type for second order loop filters.
"""
abstract type AbstractSecondOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

Second order bilinear loop filter.

Uses bilinear transformation for improved frequency response.
The natural frequency is `ω₀ = bandwidth / 0.53` (Kaplan & Hegarty, Table 5.6).

$(TYPEDFIELDS)

# Example
```julia
lf = SecondOrderBilinearLF()
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct SecondOrderBilinearLF{T} <: AbstractSecondOrderLF
    "Frequency state estimate"
    x::T
end

"""
$(TYPEDEF)

Second order boxcar loop filter.

Uses boxcar (rectangular) integration for simpler implementation.
The natural frequency is `ω₀ = bandwidth / 0.53` (Kaplan & Hegarty, Table 5.6).

$(TYPEDFIELDS)

# Example
```julia
lf = SecondOrderBoxcarLF()
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct SecondOrderBoxcarLF{T} <: AbstractSecondOrderLF
    "Frequency state estimate"
    x::T
end

"""
    SecondOrderBilinearLF()

Construct a second order bilinear loop filter with zero initial state.
"""
function SecondOrderBilinearLF()
    SecondOrderBilinearLF(0.0Hz)
end

"""
    SecondOrderBoxcarLF()

Construct a second order boxcar loop filter with zero initial state.
"""
function SecondOrderBoxcarLF()
    SecondOrderBoxcarLF(0.0Hz)
end

"""
$(SIGNATURES)

Natural frequency `ω₀ = bandwidth / 0.53` of a second order loop filter with the
noise bandwidth `bandwidth` (Kaplan & Hegarty, Table 5.6).
"""
natural_frequency(state::AbstractSecondOrderLF, bandwidth) = bandwidth / 0.53

"""
$(SIGNATURES)

Propagate the second order loop filter state.

Updates the frequency state estimate using `x_next = x + Δt * ω₀² * δθ`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
New loop filter state with updated frequency estimate.
"""
function propagate_state(state::T, δθ, Δt, ω₀) where T <: AbstractSecondOrderLF
    T(state.x + Δt * ω₀^2 * δθ)
end

"""
$(SIGNATURES)

Calculate the filtered output for the second order bilinear loop filter.

Computes `x + (√2 * ω₀ + ω₀² * Δt / 2) * δθ`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
Filtered frequency estimate.
"""
function filtered_output(state::SecondOrderBilinearLF, δθ, Δt, ω₀)
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt / 2) * δθ
end

"""
$(SIGNATURES)

Calculate the filtered output for the second order boxcar loop filter.

Computes `x + (√2 * ω₀ + ω₀² * Δt) * δθ`, which includes the current `δθ` in the
integrated frequency state (Kaplan & Hegarty, Figure 5.19(b)).
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
function filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, ω₀)
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt) * δθ
end
