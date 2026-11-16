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
lf = SecondOrderBilinearLF(; frequency = 1000.0Hz)
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
lf = SecondOrderBoxcarLF(; frequency = 1000.0Hz)
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct SecondOrderBoxcarLF{T} <: AbstractSecondOrderLF
    "Frequency state estimate"
    x::T
end

"""
    SecondOrderBilinearLF(; frequency = 0.0Hz)

A second order bilinear loop filter at the initial `frequency`.
"""
SecondOrderBilinearLF(; frequency = 0.0Hz) = SecondOrderBilinearLF(_as_frequency_state(frequency))

"""
    SecondOrderBoxcarLF(; frequency = 0.0Hz)

A second order boxcar loop filter at the initial `frequency`.
"""
SecondOrderBoxcarLF(; frequency = 0.0Hz) = SecondOrderBoxcarLF(_as_frequency_state(frequency))

natural_frequency(state::AbstractSecondOrderLF, bandwidth) = bandwidth / 0.53

function propagate_state(state::T, δθ, Δt, ω₀) where T <: AbstractSecondOrderLF
    T(state.x + Δt * ω₀^2 * δθ)
end

function filtered_output(state::SecondOrderBilinearLF, δθ, Δt, ω₀)
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt / 2) * δθ
end

# The boxcar integrator includes the current `δθ` (Kaplan & Hegarty, Figure 5.19(b)).
function filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, ω₀)
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt) * δθ
end
