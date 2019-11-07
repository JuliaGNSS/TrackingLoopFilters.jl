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

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt`
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractSecondOrderLF
    ω₀ = bandwidth * 1.89
    T(state.x + Δt * ω₀^2 * δθ)
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt`
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
"""
function get_filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + (sqrt(2) * ω₀ + ω₀^2 * Δt / 2) * δθ
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval in seconds `Δt`
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
"""
function get_filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = bandwidth * 1.89
    state.x + sqrt(2) * ω₀ * δθ
end
