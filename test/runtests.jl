using Aqua
using LinearAlgebra
using Test
using TrackingLoopFilters

import Unitful: Hz, s

@testset "Aqua" begin
    Aqua.test_all(TrackingLoopFilters)
end

@testset "First Order Loop Filter" begin
    bandwidth = 1Hz
    loop_filter = @inferred FirstOrderLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out == 4Hz

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out == 8Hz

    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out == 12Hz
end

@testset "Second Order Boxcar Loop Filter" begin
    bandwidth = 2Hz * 0.53
    loop_filter = @inferred SecondOrderBoxcarLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out ≈ 2Hz * sqrt(2)

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ 8.0Hz + 4Hz * sqrt(2)


    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ 24.0Hz + 6Hz * sqrt(2)
end

@testset "Second Order Bilinear Loop Filter" begin
    bandwidth = 2Hz * 0.53
    loop_filter = @inferred SecondOrderBilinearLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out ≈ 4Hz + 2Hz * sqrt(2)

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ 16.0Hz + 4Hz * sqrt(2) #21.656


    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ 36.0Hz + 6Hz * sqrt(2)
end

@testset "Third Order Boxcar Loop Filter" begin
    bandwidth = 2Hz * 0.7845
    loop_filter = @inferred ThirdOrderBoxcarLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out ≈ 4.8Hz

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ 8.8Hz + 2Hz *4.8

    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ 58.4Hz + 3Hz * 4.8
end

@testset "Third Order Bilinear Loop Filter" begin
    bandwidth = 2Hz * 0.7845
    loop_filter = @inferred ThirdOrderBilinearLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out ≈ 17.2Hz

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ (24.8 + 2 * 17.2) * Hz

    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ (106.4 + 3 * 17.2) * Hz
end

@testset "Third Order Assisted Bilinear Loop Filter" begin
    bandwidth = (2Hz * 0.7845, 0.5Hz * 0.53)
    loop_filter = @inferred ThirdOrderAssistedBilinearLF()
    
    out = get_filtered_output(loop_filter, [1.0 1.0Hz], 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, [1.0 1.0Hz], 2s, bandwidth)
    @test out ≈ (17.2 + 0.9571067811865461)Hz

    out = get_filtered_output(loop_filter, [2.0 2.0Hz], 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, [2.0 2.0Hz], 2s, bandwidth)
    @test out ≈ 63.02842712474619Hz

    out = get_filtered_output(loop_filter, [3.0 3.0Hz], 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, [3.0 3.0Hz], 2s, bandwidth)
    @test out ≈ 167.61396103067892Hz
end

@testset "Third Order Assisted Bilinear Loop Filter with separate bandwidths" begin
    Δt = 0.001s
    δθ = (0.1, 2.0Hz)
    # A single bandwidth is no longer supported.
    @test_throws MethodError propagate(ThirdOrderAssistedBilinearLF(), δθ, Δt, 18Hz)
    @test_throws MethodError get_filtered_output(ThirdOrderAssistedBilinearLF(), δθ, Δt, 18Hz)

    # A zero low bandwidth is the unassisted third order filter.
    assisted, plain = ThirdOrderAssistedBilinearLF(), ThirdOrderBilinearLF()
    for k in 1:5
        out_a, assisted = @inferred filter_loop(assisted, (0.01k, 3.0Hz), Δt, (10Hz, 0.0Hz))
        out_p, plain = filter_loop(plain, 0.01k, Δt, 10Hz)
        @test out_a ≈ out_p
        @test assisted.x1 ≈ plain.x1 && assisted.x2 ≈ plain.x2
    end

    # The low bandwidth alone sets the frequency path: ω₀_assist = bandwidth_low / 0.53.
    lf = propagate(ThirdOrderAssistedBilinearLF(), (0.0, 1.0Hz), 1s, (0.0Hz, 5.3Hz))
    @test lf.x1 ≈ sqrt(2) * 10Hz * 1.0Hz * 1s
    @test lf.x2 ≈ (10Hz)^2 * 1.0Hz * 1s
end

@testset "Filter" begin
    bandwidth = 1Hz
    loop_filter = @inferred FirstOrderLF()

    out, loop_filter = @inferred filter_loop(loop_filter, 1.0, 2s, bandwidth)
    @test out == 4Hz

    out, loop_filter = @inferred filter_loop(loop_filter, 2.0, 2s, bandwidth)
    @test out == 8Hz
end
