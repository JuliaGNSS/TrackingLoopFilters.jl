
abstract type   AbstractLoopFilter end

mutable struct FirstOrderLF <: AbstractLoopFilter
    x::Float64   
end

mutable struct SecondOrderBilinearLF <: AbstractLoopFilter
    x::Float64
end

mutable struct SecondOrderBoxcarLF <: AbstractLoopFilter
    x::Float64
end


mutable struct ThirdOrderBilinearLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end

mutable struct ThirdOrderBoxcarLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end