
using LinearAlgebra
using Test
using TrackingLoopFilters


import Unitful: MHz, kHz, Hz, s, ms


@testset "First Order Loop Filter" begin
    bandwidth   =   1Hz
    FrstOrdLF   =   FirstOrderLF(0)

    out = get_filtered_output(FrstOrdLF, 1.0, 2s, bandwidth)
    FrstOrdLF = propagate(FrstOrdLF, 1.0, 2s, bandwidth)
    @test out == 4Hz

    out = get_filtered_output(FrstOrdLF, 2.0, 2s, bandwidth)
    FrstOrdLF = propagate(FrstOrdLF, 2.0, 2s, bandwidth)
    @test out == 8Hz

    out = get_filtered_output(FrstOrdLF, 3.0, 2s, bandwidth)
    FrstOrdLF = propagate(FrstOrdLF, 3.0, 2s, bandwidth)
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
    @test out == 4Hz + 2 * sqrt(2)Hz

    out = get_filtered_output(SecOrdBiLF, 2.0, 2s, bandwidth)
    SecOrdBiLF = propagate(SecOrdBiLF, 2.0, 2s, bandwidth)
    @test out == 16.0Hz + 4Hz * sqrt(2) #21.656


    out = get_filtered_output(SecOrdBiLF, 3.0, 2s, bandwidth)
    SecOrdBiLF = propagate(SecOrdBiLF, 3.0, 2s, bandwidth)
    @test out == 36.0Hz + 6Hz * sqrt(2)
end


@testset "Third Order Boxcar Loop Filter" begin
    bandwidth = 2Hz / 1.2
    ThrdOrdBcLF = ThirdOrderBoxcarLF([0 ; 0])

    out = get_filtered_output(ThrdOrdBcLF, 1.0, 2s, bandwidth)
    ThrdOrdBcLF = propagate(ThrdOrdBcLF, 1.0, 2s, bandwidth)
    @test out == 4.8Hz

    out = get_filtered_output(ThrdOrdBcLF, 2.0, 2s, bandwidth)
    ThrdOrdBcLF = propagate(ThrdOrdBcLF, 2.0, 2s, bandwidth)
    @test out == 8.8Hz + 2Hz *4.8

    out = get_filtered_output(ThrdOrdBcLF, 3.0, 2s, bandwidth)
    ThrdOrdBcLF = propagate(ThrdOrdBcLF, 3.0, 2s, bandwidth)
    @test out == 58.4Hz + 3Hz * 4.8
end


@testset "Third Order Bilinear Loop Filter" begin
    bandwidth = 2Hz / 1.2
    ThrdOrdBiLF = ThirdOrderBilinearLF([0 ; 0])

    out = get_filtered_output(ThrdOrdBiLF, 1.0, 2s, bandwidth)
    ThrdOrdBiLF = propagate(ThrdOrdBiLF, 1.0, 2s, bandwidth)
    @test out == 17.2Hz

    out = get_filtered_output(ThrdOrdBiLF, 2.0, 2s, bandwidth)
    ThrdOrdBiLF = propagate(ThrdOrdBiLF, 2.0, 2s, bandwidth)
    @test out == (24.8 + 2 * 17.2) * Hz

    out = get_filtered_output(ThrdOrdBiLF, 3.0, 2s, bandwidth)
    ThrdOrdBiLF = propagate(ThrdOrdBiLF, 3.0, 2s, bandwidth)
    @test out ==  (106.4 + 3 * 17.2) * Hz
end
