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
    @test out ≈ 8.0Hz + 2Hz * sqrt(2)

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ 24.0Hz + 4Hz * sqrt(2)


    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ 48.0Hz + 6Hz * sqrt(2)
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
    @test out ≈ 45.6Hz

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ (40.8 + 32 + 2 * 45.6) * Hz

    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ (154.4 + 96 + 3 * 45.6) * Hz
end

@testset "Third Order Bilinear Loop Filter" begin
    bandwidth = 2Hz * 0.7845
    loop_filter = @inferred ThirdOrderBilinearLF()

    out = @inferred get_filtered_output(loop_filter, 1.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 1.0, 2s, bandwidth)
    @test out ≈ 17.2Hz

    out = @inferred get_filtered_output(loop_filter, 2.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 2.0, 2s, bandwidth)
    @test out ≈ (24.8 + 16 + 2 * 17.2) * Hz

    out = @inferred get_filtered_output(loop_filter, 3.0, 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, 3.0, 2s, bandwidth)
    @test out ≈ (106.4 + 48 + 3 * 17.2) * Hz
end

@testset "Third Order Assisted Bilinear Loop Filter" begin
    bandwidth = (2Hz * 0.7845, 0.5Hz * 0.53)
    loop_filter = @inferred ThirdOrderAssistedBilinearLF()

    out = @inferred get_filtered_output(loop_filter, (1.0, 1.0Hz), 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, (1.0, 1.0Hz), 2s, bandwidth)
    @test out ≈ (17.2 + 0.9571067811865461)Hz

    out = @inferred get_filtered_output(loop_filter, (2.0, 2.0Hz), 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, (2.0, 2.0Hz), 2s, bandwidth)
    @test out ≈ 79.52842712474619Hz

    out = @inferred get_filtered_output(loop_filter, (3.0, 3.0Hz), 2s, bandwidth)
    loop_filter = @inferred propagate(loop_filter, (3.0, 3.0Hz), 2s, bandwidth)
    @test out ≈ 217.1139610306789Hz
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
    @test lf.x1 ≈ (sqrt(2) * 10Hz + (10Hz)^2 * 1s / 2) * 1.0Hz * 1s
    @test lf.x2 ≈ (10Hz)^2 * 1.0Hz * 1s
end

# Integrators of Kaplan & Hegarty, Figure 5.19: the boxcar integrator outputs the
# new accumulator, the bilinear one the mean of the old and the new.
boxcar(A, A_next) = A_next
bilinear(A, A_next) = (A + A_next) / 2

# Second order filter of Kaplan & Hegarty, Figure 5.20(b), written out with the
# accumulator of its integrator.
function kaplan_second_order(δθs, Δt, ω₀, integrator)
    A = 0.0Hz
    map(δθs) do δθ
        A_next = A + Δt * ω₀^2 * δθ
        out = integrator(A, A_next) + sqrt(2) * ω₀ * δθ
        A = A_next
        out
    end
end

# Third order filter of Kaplan & Hegarty, Figure 5.20(c), written out with the
# accumulators of its two integrators. `δf` and `ω₀f` feed the frequency error
# into the FLL-assisted variant of Figure 5.21(b).
function kaplan_third_order(δθs, δfs, Δt, ω₀, ω₀f, integrator)
    A1, A2 = 0.0Hz^2, 0.0Hz
    map(zip(δθs, δfs)) do (δθ, δf)
        A1_next = A1 + Δt * (ω₀^3 * δθ + ω₀f^2 * δf)
        y1 = integrator(A1, A1_next) + 1.1 * ω₀^2 * δθ + sqrt(2) * ω₀f * δf
        A2_next = A2 + Δt * y1
        out = integrator(A2, A2_next) + 2.4 * ω₀ * δθ
        A1, A2 = A1_next, A2_next
        out
    end
end

@testset "Loop Filters follow Kaplan & Hegarty" begin
    Δt = 0.004s
    δθs = [0.1, -0.05, 0.2, 0.0, 0.15, -0.1]
    δfs = [2.0, 1.0, -1.5, 0.5, 0.0, 3.0] .* Hz

    filters = (
        (SecondOrderBoxcarLF(), δθs, 4Hz, kaplan_second_order(δθs, Δt, 4Hz / 0.53, boxcar)),
        (SecondOrderBilinearLF(), δθs, 4Hz, kaplan_second_order(δθs, Δt, 4Hz / 0.53, bilinear)),
        (ThirdOrderBoxcarLF(), δθs, 15Hz, kaplan_third_order(δθs, zero(δfs), Δt, 15Hz / 0.7845, 0.0Hz, boxcar)),
        (ThirdOrderBilinearLF(), δθs, 15Hz, kaplan_third_order(δθs, zero(δfs), Δt, 15Hz / 0.7845, 0.0Hz, bilinear)),
        (ThirdOrderAssistedBilinearLF(), collect(zip(δθs, δfs)), (15Hz, 4Hz),
            kaplan_third_order(δθs, δfs, Δt, 15Hz / 0.7845, 4Hz / 0.53, bilinear)),
    )
    for (lf, inputs, bandwidth, expected) in filters
        for (δθ, out_expected) in zip(inputs, expected)
            out, lf = @inferred filter_loop(lf, δθ, Δt, bandwidth)
            @test out ≈ out_expected
        end
    end
end

@testset "Filter" begin
    bandwidth = 1Hz
    loop_filter = @inferred FirstOrderLF()

    out, loop_filter = @inferred filter_loop(loop_filter, 1.0, 2s, bandwidth)
    @test out == 4Hz

    out, loop_filter = @inferred filter_loop(loop_filter, 2.0, 2s, bandwidth)
    @test out == 8Hz
end
