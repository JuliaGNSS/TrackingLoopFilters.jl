"""
$(TYPEDEF)

Abstract base type for first order loop filters.
"""
abstract type AbstractFirstOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

First order loop filter (proportional only).

A stateless filter that provides proportional gain without integration.
The natural frequency scaling factor is 4.0.

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

Propagate the first order loop filter state.

Since the first order filter is stateless, this returns the input state unchanged.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
The unchanged loop filter state.
"""
function propagate(state::FirstOrderLF, δθ, Δt, bandwidth)
    state
end

"""
$(SIGNATURES)

Calculate the filtered output for the first order loop filter.

Computes `ω₀ * δθ` where `ω₀ = bandwidth * 4.0`.

# Arguments
- `state`: Current loop filter state
- `δθ`: Phase discriminator output
- `Δt`: Integration time
- `bandwidth`: Loop bandwidth

# Returns
Filtered frequency estimate.
"""
function get_filtered_output(state::FirstOrderLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 4.0
    ω₀ * δθ
end
