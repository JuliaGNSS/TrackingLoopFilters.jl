abstract type AbstractFirstOrderLF <: AbstractLoopFilter end

struct FirstOrderLF <: AbstractFirstOrderLF
end

"""
$(SIGNATURES)

Propagates the state of the loop filter.
"""
function propagate(state::FirstOrderLF, δθ, Δt, bandwidth)
    state
end

"""
$(SIGNATURES)

Calculates the output of the loop filter.
"""
function get_filtered_output(state::FirstOrderLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 4.0
    ω₀ * δθ
end
