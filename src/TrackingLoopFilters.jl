module TrackingLoopFilters

    using DocStringExtensions
    using LinearAlgebra
    using StaticArrays

    import Unitful: MHz, kHz, Hz, s, ms, upreferred



    export
        propagate,
        get_filtered_output,
        FirstOrderLF,
        SecondOrderBilinearLF,
        SecondOrderBoxcarLF,
        ThirdOrderBilinearLF,
        ThirdOrderBoxcarLF


    abstract type AbstractLoopFilter end

    include("FirstOrderLF.jl")
    include("SecondOrderLF.jl")
    include("ThirdOrderLF.jl") 

    #function loop_filter(state::AbstractLoopFilter, δθ, Δt, bandwidth)
    #    Δt_sec = Float64(upreferred(Δt/s))
    #    out = get_filtered_output(state, δθ, Δt_sec, bandwidth) * Hz 
    #    state.x = propagate(state, δθ, Δt_sec, bandwidth).x
    #    out
    #end
end








 




    






