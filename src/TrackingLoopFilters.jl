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

end








 




    






