mutable struct ThirdOrderBilinearLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end

mutable struct ThirdOrderBoxcarLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end


"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.2
    F = @SMatrix [1.0 Δt; 0.0 1.0]
    L = @SVector [Δt * 1.1 * ω₀^2, Δt * ω₀^3]
    ThirdOrderBoxcarLF(F * state.x + L * δθ)
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.2
    F =  @SMatrix [1.0 Δt; 0.0 1.0]
    L =  @SVector [Δt * 1.1 * ω₀^2, Δt * ω₀^3]
    ThirdOrderBilinearLF(F * state.x + L * δθ)
end





"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
"""    
function get_filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    C = @SVector [1.0, Δt / 2]
    D = 2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4
    dot(C , state.x) + D * δθ
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    C = @SVector [1.0, 0.0]
    D = 2.4 * ω₀
    dot(C , state.x) + D * δθ
end  


