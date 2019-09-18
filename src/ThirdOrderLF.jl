using StaticArrays, LinearAlgebra

include("FilterStructs.jl")

function propagate(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.2
    F = @SMatrix [1.0 Δt; 0.0 1.0]
    L = @SVector [Δt * 1.1 * ω₀^2, Δt * ω₀^3]

    println(F, state.x, L, δθ)
    ThirdOrderBoxcarLF(F * state.x + L * δθ)
end

function propagate(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.2
    F =  @SMatrix [1.0 Δt; 0.0 1.0]
    L =  @SVector [Δt * 1.1 * ω₀^2, Δt * ω₀^3]
    ThirdOrderBilinearLF(F * state.x + L * δθ)
end




    
function filtered_output(state::ThirdOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    C = @SVector [1.0, Δt / 2]
    D = 2.4 * ω₀ + 1.1 * ω₀^2 * Δt / 2 + ω₀^3 * Δt^2 / 4
    dot(C , state.x) + D * δθ
end

function filtered_output(state::ThirdOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.2
    C = @SVector [1.0, 0.0]
    D = 2.4 * ω₀
    dot(C , state.x) + D * δθ
end  


