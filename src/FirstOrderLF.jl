"""
$(TYPEDEF)

Abstract base type for first order loop filters.
"""
abstract type AbstractFirstOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

First order loop filter (proportional only).

A stateless filter that provides proportional gain without integration.
The natural frequency is `ω₀ = bandwidth / 0.25` (Kaplan & Hegarty, Table 5.6).

# Example
```julia
lf = FirstOrderLF()
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct FirstOrderLF <: AbstractFirstOrderLF
end

"""
$(SIGNATURES)

Natural frequency `ω₀ = bandwidth / 0.25` of a first order loop filter with the
noise bandwidth `bandwidth` (Kaplan & Hegarty, Table 5.6).
"""
natural_frequency(state::AbstractFirstOrderLF, bandwidth) = bandwidth / 0.25

"""
$(SIGNATURES)

Propagate the first order loop filter state.

Since the first order filter is stateless, this returns the input state unchanged.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
The unchanged loop filter state.
"""
function propagate_state(state::FirstOrderLF, δθ, Δt, ω₀)
    state
end

"""
$(SIGNATURES)

Calculate the filtered output for the first order loop filter.

Computes `ω₀ * δθ`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `ω₀`: Natural frequency

# Returns
Filtered frequency estimate.
"""
function filtered_output(state::FirstOrderLF, δθ, Δt, ω₀)
    ω₀ * δθ
end
