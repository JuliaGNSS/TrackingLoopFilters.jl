
abstract type AbstractThirdOrderLF <: AbstractLoopFilter end
abstract type AbstractThirdOrderAssistedLF <: AbstractLoopFilter end

struct ThirdOrderBilinearLF{T1,T2} <: AbstractThirdOrderLF
    x1::T1
    x2::T2
end

struct ThirdOrderAssistedBilinearLF{T1,T2} <: AbstractThirdOrderAssistedLF
    x1::T1
    x2::T2
end

struct ThirdOrderBoxcarLF{T1,T2} <: AbstractThirdOrderLF
    x1::T1
    x2::T2
end

function ThirdOrderBilinearLF()
    ThirdOrderBilinearLF(0.0Hz, 0.0Hz^2)
end

function ThirdOrderAssistedBilinearLF()
    ThirdOrderAssistedBilinearLF(0.0Hz, 0.0Hz^2)
end

function ThirdOrderBoxcarLF()
    ThirdOrderBoxcarLF(0.0Hz, 0.0Hz^2)
end

"""
$(SIGNATURES)

Propagates the state of the loop filter.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractThirdOrderLF
    ω₀ = bandwidth * 1.2
    T(state.x1 + Δt * state.x2 + 1.1 * Δt * ω₀^2 * δθ, state.x2 + Δt * ω₀^3 * δθ)
end

"""
$(SIGNATURES)

Propagates the state of the assisted loop filter.
"""
function propagate(state::T, δθ, Δt, bandwidth) where T <: AbstractThirdOrderAssistedLF
    ω₀ = bandwidth * 1.2
    ω₀_assist = ω₀ / 4
    T(state.x1 + Δt * state.x2 + 1.1 * Δt * ω₀^2 * δθ[1] + Δt * sqrt(2) * ω₀_assist * δθ[2], state.x2 + Δt * ω₀^3 * δθ[1] + Δt * ω₀_assist^2 * δθ[2])
end

"""
$(SIGNATURES)

Calculates the output of the loop filter.
"""
function get_filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀= bandwidth * 1.2
    state.x1 + Δt / 2 * state.x2 + (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ
end

"""
$(SIGNATURES)

Calculates the output of the third order loop filter assisted by a second order loop filter.
"""
function get_filtered_output(state::ThirdOrderAssistedBilinearLF, δθ, Δt, bandwidth)
    ω₀= bandwidth * 1.2
    ω₀_assist = ω₀ / 4
    state.x1 + Δt / 2 * state.x2 +
        (2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4) * δθ[1] +
        (sqrt(2) * ω₀_assist * Δt / 2 + ω₀_assist^2 * Δt^2 / 4) * δθ[2]
end

"""
$(SIGNATURES)

Calculates the output of the loop filter.
"""
function get_filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= bandwidth * 1.2
    state.x1 + 2.4 * ω₀ * δθ
end
