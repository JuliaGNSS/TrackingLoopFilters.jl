
    using StaticArrays, LinearAlgebra
    import Unitful: MHz, kHz, Hz, s, ms, upreferred


    abstract type AbstractLoopFilter end



    """

    """
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




    """
    #TODO:
    # - Wo einbinden?
    # - genaue Funktion?
    """

    """
    #TODO: 
    # [x] Für andere Structs noch fortführen
    # - funktion einbinden
    """


     function loop_filter(state::T, δθ, Δt, bandwidth) where T<:AbstractLoopFilter
         ω₀ = Float64(bandwidth / Hz) * 4.0
         Δt_sec = Float64(upreferred(Δt/s))
         propagate(state, δθ, Δt_sec, bandwidth)
         out = filtered_output(state, δθ, Δt_sec, bandwidth) * Hz


     end


    # function loop_filter(state::SecondOrderBilinearLF, Δt, δθ, bandwidth)
    #     ω₀ = Float64(bandwidth / Hz) * 1.89
    #     propagate(state, δθ, Δt, bandwidth)
    #     filtered_output(state, δθ, Δt, bandwidth)
    # end

    # function loop_filter(state::SecondOrderBoxcarLF, Δt, δθ, bandwidth)
    #     ω₀ = Float64(bandwidth / Hz) * 1.89
    #     next_x= propagate(state, δθ, Δt, bandwidth).x
    #     filtered_output(state, δθ, Δt, bandwidth)
    # end

    # function loop_filter(state::ThirdOrderBilinearLF, Δt, δθ, bandwidth)
    #     ω₀ = Float64(bandwidth / Hz) * 1.2
    #     propagate(state, δθ, Δt, bandwidth)
    #     filtered_output(state, δθ, Δt, bandwidth)
    # end

    # function loop_filter(state::ThirdOrderBoxcarLF, Δt, δθ, bandwidth)    
    #     ω₀ = Float64(bandwidth / Hz) * 1.2
    #     propagate(state, δθ, Δt, bandwidth)
    #     filtered_output(state, δθ, Δt, bandwidth)
    # end









    """
    # i.O. wenn ich zusätzliche Übergabeparameter für Firstorder festlege? -> nur eine Loop Filter Funktion
    """

    function propagate(state::FirstOrderLF, δθ, Δt, bandwidth)
        state
    end

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





    """

    """
    function filtered_output(state::FirstOrderLF, δθ, Δt, bandwidth)
        ω₀ = Float64(bandwidth/Hz) * 4.0
        ω₀ * δθ
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


    






