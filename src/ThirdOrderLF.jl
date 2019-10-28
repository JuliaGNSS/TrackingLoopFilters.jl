
abstract type AbstractThirdOrderLF <: AbstractLoopFilter end

struct ThirdOrderBilinearLF <: AbstractThirdOrderLF
    x::SVector{2, Float64}
end

struct ThirdOrderBoxcarLF <: AbstractThirdOrderLF
    x::SVector{2, Float64}
end

function ThirdOrderBilinearLF()
    ThirdOrderBilinearLF(@SVector[0 , 0] )
end

function ThirdOrderBoxcarLF()
    ThirdOrderBoxcarLF(@SVector[0 , 0])
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::T, δθ, Δt, bandwidth) where T<:AbstractThirdOrderLF
    ω₀ = Float64(bandwidth/Hz) * 1.2
    Δt = Δt/s
    F = @SMatrix [1.0 Δt; 0.0 1.0]
    L = @SVector [Δt * 1.1 * ω₀^2, Δt * ω₀^3]
    T(F * state.x + L * δθ)
end





"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
"""    
function get_filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    Δt = Δt/s
    C = @SVector [1.0, Δt / 2]
    D = 2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4
    (dot(C , state.x) + D * δθ) * Hz
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    Δt = Δt/s
    C = @SVector [1.0, 0.0]
    D = 2.4 * ω₀
    (dot(C , state.x) + D * δθ) * Hz
end  


