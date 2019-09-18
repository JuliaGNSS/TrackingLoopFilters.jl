import Unitful: MHz, kHz, Hz, s, ms
include("FilterStructs.jl")

function propagate(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    F = 1.0
    L = Δt * ω₀^2
    SecondOrderBoxcarLF(F * state.x + L * δθ)
end

function propagate(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 1.89
    F = 1.0
    L = Δt * ω₀^2
    SecondOrderBilinearLF(F * state.x + L * δθ)
end





function filtered_output(state::SecondOrderBilinearLF, δθ, Δt, bandwidth)
    ω₀=Float64(bandwidth/Hz) * 1.89
    C = 1.0
    D = sqrt(2) * ω₀ + ω₀^2 * Δt / 2
    C * state.x + D * δθ
end

function filtered_output(state::SecondOrderBoxcarLF, δθ, Δt, bandwidth)
    ω₀= Float64(bandwidth/Hz) * 1.89
    C = 1.0 
    D = sqrt(2) * ω₀
    C * state.x + D * δθ
end

