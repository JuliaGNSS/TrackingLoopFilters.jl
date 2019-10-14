
abstract type AbstractSecondOrderLF <: AbstractLoopFilter end

struct SecondOrderBilinearLF <: AbstractSecondOrderLF
    x::Float64
end

struct SecondOrderBoxcarLF <: AbstractSecondOrderLF
    x::Float64
end


"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::AbstractSecondOrderLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    Δt = Δt/s
    F = 1.0
    L = Δt * ω₀^2

    if typeof(state) == SecondOrderBoxcarLF
        return SecondOrderBoxcarLF(F * state.x + L * δθ)

    elseif typeof(state) == SecondOrderBilinearLF
        return SecondOrderBilinearLF(F * state.x + L * δθ)

    end
end





"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    Δt = Δt/s
    C = 1.0
    D = sqrt(2) * ω₀ + ω₀^2 * Δt / 2
    (C * state.x + D * δθ)*Hz
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval in seconds `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    Δt = Δt/s
    C = 1.0 
    D = sqrt(2) * ω₀
    (C * state.x + D * δθ)*Hz
end

