
include("FilterStructs.jl")

function propagate(state::FirstOrderLF, δθ, Δt, bandwidth)
    state
end


function filtered_output(state::FirstOrderLF, δθ, Δt, bandwidth)
    ω₀ = Float64(bandwidth/Hz) * 4.0
    ω₀ * δθ
end

