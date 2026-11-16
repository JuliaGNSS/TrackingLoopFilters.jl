"""
$(TYPEDEF)

Abstract base type for first order loop filters.
"""
abstract type AbstractFirstOrderLF <: AbstractLoopFilter end

"""
$(TYPEDEF)

First order loop filter (proportional only).

Provides proportional gain without integration, around a constant frequency
that propagation leaves unchanged.
The natural frequency is `ω₀ = bandwidth / 0.25` (Kaplan & Hegarty, Table 5.6).

$(TYPEDFIELDS)

# Example
```julia
lf = FirstOrderLF()
lf = FirstOrderLF(; frequency = 1000.0Hz)
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```
"""
struct FirstOrderLF{T} <: AbstractFirstOrderLF
    "Constant frequency"
    x::T
end

"""
    FirstOrderLF(; frequency = 0.0Hz)

A first order loop filter around the constant `frequency`.
"""
FirstOrderLF(; frequency = 0.0Hz) = FirstOrderLF(_as_frequency_state(frequency))

natural_frequency(state::AbstractFirstOrderLF, bandwidth) = bandwidth / 0.25

# The frequency is constant.
function propagate_state(state::FirstOrderLF, δθ, Δt, ω₀)
    state
end

function filtered_output(state::FirstOrderLF, δθ, Δt, ω₀)
    state.x + ω₀ * δθ
end

