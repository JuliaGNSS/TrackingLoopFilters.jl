abstract type AbstractSecondOrderLF <: AbstractLoopFilter end

struct SecondOrderBilinearLF{T} <: AbstractSecondOrderLF
    x::T
end

struct SecondOrderBoxcarLF{T} <: AbstractSecondOrderLF
    x::T
end

function SecondOrderBilinearLF()
    SecondOrderBilinearLF(0.0Hz)
end

function SecondOrderBoxcarLF()
    SecondOrderBoxcarLF(0.0Hz)
end

"""
$(SIGNATURES)

Propagates the state of the loop filter.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractSecondOrderLF
    ω₀ = bandwidth * 1.89
    T(state.x + Δt * ω₀^2 * δθ)
end

"""
$(SIGNATURES)

Calculates the output of the loop filter.
"""
function get_filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt / 2) * δθ
end

"""
$(SIGNATURES)

Calculates the output of the loop filter.
"""
function get_filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + sqrt(2) * ω₀ * δθ
end
