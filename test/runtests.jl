
using LinearAlgebra
using Test
using TrackingLoopFilters


import Unitful: MHz, kHz, Hz, s, ms


@testset "First Order Loop Filter" begin
    bandwidth   =   1Hz
    FrstOrdLF   =   FirstOrderLF(0)
    out = loop_filter(FrstOrdLF, 1.0, 2s, bandwidth)
    @test out == 4Hz
    out = loop_filter(FrstOrdLF, 2.0, 2s, bandwidth)
    @test out == 8Hz
    out = loop_filter(FrstOrdLF, 3.0, 2s, bandwidth)
    @test out == 12Hz
end


@testset "Second Order Boxcar Loop Filter" begin
    bandwidth = 2Hz / 1.89
    SecOrdBcLF = SecondOrderBoxcarLF(0)

    out = get_filtered_output(SecOrdBcLF, 1.0, 2s, bandwidth)
    SecOrdBcLF = propagate(SecOrdBcLF, 1.0, 2s, bandwidth)
    @test out == 2 * sqrt(2)Hz

    out = get_filtered_output(SecOrdBcLF, 2.0, 2s, bandwidth)
    SecOrdBcLF = propagate(SecOrdBcLF, 2.0, 2s, bandwidth)
    @test out == 8.0Hz + 4Hz * sqrt(2)


    out = get_filtered_output(SecOrdBcLF, 3.0, 2s, bandwidth)
    SecOrdBcLF = propagate(SecOrdBcLF, 3.0, 2s, bandwidth)
    @test out == 24.0Hz + 6Hz * sqrt(2)
end

@testset "Second Order Bilinear Loop Filter" begin
    bandwidth = 2Hz / 1.89
    SecOrdBiLF = SecondOrderBilinearLF(0)

    out = get_filtered_output(SecOrdBiLF, 1.0, 2s, bandwidth)
    SecOrdBiLF = propagate(SecOrdBiLF, 1.0, 2s, bandwidth)
    @test out == 2 * sqrt(2)Hz

    out = get_filtered_output(SecOrdBiLF, 2.0, 2s, bandwidth)
    SecOrdBiLF = propagate(SecOrdBiLF, 2.0, 2s, bandwidth)
    @test out == 8.0Hz + 4Hz * sqrt(2)


    out = get_filtered_output(SecOrdBiLF, 3.0, 2s, bandwidth)
    SecOrdBiLF = propagate(SecOrdBiLF, 3.0, 2s, bandwidth)
    @test out == 24.0Hz + 6Hz * sqrt(2)
end



@testset "Third Order Boxcar Loop Filter" begin
    bandwidth = 2Hz / 1.2
    ThrdOrdBcLF = ThirdOrderBoxcarLF([0 ; 0])
    out = loop_filter(ThrdOrdBcLF, 1.0, 2s, bandwidth)
    @test out == 4.8Hz
    out = loop_filter(ThrdOrdBcLF, 2.0, 2s, bandwidth)
    @test out == 8.8Hz + 2Hz *4.8
    out = loop_filter(ThrdOrdBcLF, 3.0, 2s, bandwidth)
    @test out == 58.4Hz + 3Hz * 4.8
end

