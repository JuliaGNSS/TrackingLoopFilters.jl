using StaticArrays

abstract type AbstractLoopFilter end

struct FirstOrderLF <: AbstractLoopFilter   
end

struct SecondOrderBilinearLF <: AbstractLoopFilter
    x::Float64
end

struct SecondOrderBoxcarLF <: AbstractLoopFilter
    x::Float64
end


struct ThirdOrderBilinearLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end

struct ThirdOrderBoxcarLF <: AbstractLoopFilter
    x::SVector{2, Float64}
end