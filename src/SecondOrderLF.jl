"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    F = 1.0
    L = Δt * ω₀^2
    SecondOrderBoxcarLF(F * state.x + L * δθ)
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `F` and `L` (Transition Matrix and Filter gain Matrix)
matrices to calculate the initial state vector `x` and create a new object
of the same type with new state
"""
function propagate(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    F = 1.0
    L = Δt * ω₀^2
    SecondOrderBilinearLF(F * state.x + L * δθ)
end




"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀=Float64(bandwidth/Hz) * 1.89
    C = 1.0
    D = sqrt(2) * ω₀ + ω₀^2 * Δt / 2
    C * state.x + D * δθ
end

"""
$(SIGNATURES)

Uses the current state, the discriminator output `δθ`, the loop update time interval `Δt` 
and the loop bandwidth `bandwidth` to set up the `C` and `D` (Transition Matrix and Filter gain Matrix)
matrices to calculate the the system output
""" 
function get_filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.89
    C = 1.0 
    D = sqrt(2) * ω₀
    C * state.x + D * δθ
end

