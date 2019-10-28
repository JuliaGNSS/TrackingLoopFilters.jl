abstract type AbstractFirstOrderLF <: AbstractLoopFilter end

struct FirstOrderLF <: AbstractFirstOrderLF
    
end

"""
$(SIGNATURES)

returns the state again. 
"""
function propagate(state::FirstOrderLF, δθ, Δt, bandwidth)
    state
end



"""
$(SIGNATURES)
     
Uses the discriminator output `δθ` and the loop bandwidth `bandwidth` 
to calculate dthe product of ω₀ and δθ as system output
"""
function get_filtered_output(state::FirstOrderLF, δθ, Δt, bandwidth)
   
    ω₀ = Float64(bandwidth/Hz) * 4.0
    (ω₀ * δθ) * Hz
end

