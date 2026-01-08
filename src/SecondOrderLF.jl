"""
$(TYPEDEF)

Abstract base type for second order loop filters.
"""
abstract type AbstractSecondOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

Second order bilinear loop filter.

Uses bilinear transformation for improved frequency response.
The natural frequency scaling factor is 1.89.

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
The natural frequency scaling factor is 1.89.

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

Propagate the second order loop filter state.

Updates the frequency state estimate using `x_next = x + Δt * ω₀² * δθ`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
New loop filter state with updated frequency estimate.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractSecondOrderLF
    ω₀ = bandwidth * 1.89
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
- `bandwidth`: Loop bandwidth

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt / 2) * δθ
end

"""
$(SIGNATURES)

Calculate the filtered output for the second order boxcar loop filter.

Computes `x + √2 * ω₀ * δθ`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + sqrt(2) * ω₀ * δθ
end
