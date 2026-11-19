using Aqua
using Test
using TrackingLoopFilters
using TrackingLoopFilters: get_filtered_output, propagate

import Unitful: DimensionError, Hz, s, unit, @u_str

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
    # The assisted filter takes a Tuple bandwidth; a single one raises a MethodError.
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

@testset "Propagating keeps the state types" begin
    for (lf, δθ, bandwidth) in (
        (ThirdOrderBilinearLF(), 0.1, 15Hz),
        (ThirdOrderBoxcarLF(), 0.1, 15Hz),
        (ThirdOrderAssistedBilinearLF(), (0.1, 1.0Hz), (15Hz, 4Hz)),
        (SecondOrderBilinearLF(), 0.1, 15Hz),
        (SecondOrderBoxcarLF(), 0.1, 15Hz),
    )
        @test typeof(propagate(lf, δθ, 0.004s, bandwidth)) == typeof(lf)
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

# Tracks the Doppler f(t) = f0 + a·t for N integrations of T, starting at t = 0 with the
# carrier NCO at `command` and the residual phase (cycles) at `φ0`. Each output commands
# the next integration. The discriminators are the mean residual phase over the
# integration and, for an assisted filter, the mean residual frequency (Hz). Returns the
# filter, the discriminator outputs and the commands.
function track_ramp(lf, bandwidth; f0, a, T, N, command = 0.0, φ0 = 0.0)
    φ = φ0
    δθs, commands = Float64[], Float64[]
    for k in 0:N-1
        t0 = k * T
        δf = f0 + a * (t0 + T / 2) - command
        δθ = φ + (f0 + a * t0 - command) * T / 2 + a * T^2 / 6
        φ += δf * T
        input = is_assisted(lf) ? (δθ, δf * Hz) : δθ
        out, lf = filter_loop(lf, input, T * s, bandwidth)
        command = Float64(out / Hz)
        push!(δθs, δθ)
        push!(commands, command)
    end
    lf, δθs, commands
end

@testset "Frequency estimates" begin
    f0, a, T, N = 100.0, 50.0, 1e-3, 20_000
    t_end = N * T
    truth(t) = (f0 + a * t) * Hz
    for (lf, bandwidth) in (
        (ThirdOrderBilinearLF(), 10Hz),
        (ThirdOrderBoxcarLF(), 10Hz),
        (ThirdOrderAssistedBilinearLF(), (10Hz, 2Hz)),
    )
        lf, = track_ramp(lf, bandwidth; f0, a, T, N)
        @test @inferred(get_frequency_rate(lf)) ≈ a * Hz / s
        @test unit(get_frequency_rate(lf)) == Hz / s
        for τ in (-T / 2, 0.0, T / 2, 0.1)
            @test @inferred(predict_frequency(lf, τ * s)) ≈ truth(t_end + τ)
            @test unit(predict_frequency(lf, τ * s)) == Hz
        end
        # Without a new discriminator output, the output after the next integration is
        # the frequency at the centre of the one after it.
        @test predict_frequency(lf, 3T / 2 * s) ≈
            get_filtered_output(lf, is_assisted(lf) ? (0.0, 0.0Hz) : 0.0, T * s, bandwidth)
    end
    for lf in (SecondOrderBilinearLF(), SecondOrderBoxcarLF())
        lf, = track_ramp(lf, 5Hz; f0, a = 0.0, T, N)
        @test predict_frequency(lf, 0.1s) ≈ f0 * Hz
        @test get_frequency_rate(lf) === 0.0Hz / s
    end

    # `τ` is a time for every filter, also below third order, which ignore it.
    for lf in (SecondOrderBilinearLF(), SecondOrderBoxcarLF(), ThirdOrderBilinearLF(),
        ThirdOrderBoxcarLF(), ThirdOrderAssistedBilinearLF())
        @test_throws DimensionError predict_frequency(lf, 0)
    end

    @test @inferred(is_assisted(ThirdOrderAssistedBilinearLF()))
    @test is_assisted(ThirdOrderAssistedBilinearLF)
    @test !is_assisted(ThirdOrderBilinearLF())
    @test !is_assisted(SecondOrderBoxcarLF)
    @test !is_assisted(FirstOrderLF())
end

# Tracks the Doppler f(t) = f0 + a·t over integrations of the lengths `Ts`, starting at
# t = 0 with the residual phase (cycles) at `φ0`. Each integration is commanded with the
# filter's prediction for its centre, as a tracking loop that knows the length of the
# next integration does. The discriminators are as in `track_ramp`, or zero for
# `open_loop`. Returns the filter after each integration and the command errors.
function track_ramp_predicted(lf, bandwidth, Ts; f0, a, φ0 = 0.0, open_loop = false)
    φ, t0 = φ0, 0.0
    lfs, command_errors = typeof(lf)[], Float64[]
    for T in Ts
        command = Float64(predict_frequency(lf, T / 2 * s) / Hz)
        δf = f0 + a * (t0 + T / 2) - command
        δθ = open_loop ? 0.0 : φ + (f0 + a * t0 - command) * T / 2 + a * T^2 / 6
        φ += δf * T
        input = is_assisted(lf) ? (δθ, (open_loop ? 0.0 : δf) * Hz) : δθ
        _, lf = filter_loop(lf, input, T * s, bandwidth)
        t0 += T
        push!(lfs, lf)
        push!(command_errors, -δf)
    end
    lfs, command_errors
end

@testset "Frequency estimates across a change of the integration time" begin
    f0, a = 100.0, 50.0
    T1, T2 = 1e-3, 20e-3
    # 1 ms integrations, then 20 ms after bit synchronisation, then 1 ms again.
    Ts = [fill(T1, 2000); fill(T2, 100); fill(T1, 100)]
    ends = cumsum(Ts)
    truth(t) = (f0 + a * t) * Hz
    # 2 Hz keep the loops stable at 20 ms.
    for (F, bandwidth) in ((ThirdOrderBilinearLF, 2Hz), (ThirdOrderBoxcarLF, 2Hz),
        (ThirdOrderAssistedBilinearLF, (2Hz, 1Hz)))
        # Without discriminator outputs, the filter extrapolates its estimates of the
        # ramp exactly, whatever the integration times.
        lfs, command_errors = track_ramp_predicted(F(f0 * Hz, a * Hz^2), bandwidth, Ts;
            f0, a, open_loop = true)
        @test maximum(abs, command_errors) < 1e-9
        for (lf, t_end, T) in zip(lfs, ends, Ts)
            @test predict_frequency(lf, 0.0s) ≈ truth(t_end) rtol = 1e-12
            @test lf.Δt == T * s
            # Without a new discriminator output, the output commands the frequency at
            # the centre of the integration after the next, if it has the same length.
            input = is_assisted(lf) ? (0.0, 0.0Hz) : 0.0
            @test filter_loop(lf, input, T * s, bandwidth)[1] ≈ truth(t_end + 3T / 2) rtol = 1e-12
        end

        # Tracking the ramp, the change of the integration time changes the mean
        # residual phase over an integration by a(T1² - T2²) / 12, which the loop
        # follows, but it does not shift the frequency estimate by (T2 - T1)·a or
        # (T2 - T1)·a / 2, 0.95 Hz or 0.475 Hz. Without moving the frequency state to
        # the new integration time, the errors below are 0.59 to 1.2 Hz after the first
        # 20 ms integration and up to 5.6 Hz later.
        lfs, command_errors = track_ramp_predicted(F(f0 * Hz, a * Hz^2), bandwidth, Ts;
            f0, a, φ0 = a * T1^2 / 12)
        @test maximum(abs, command_errors[1:2000]) < 1e-8
        switch = 2001
        @test abs(predict_frequency(lfs[switch], 0.0s) - truth(ends[switch])) < 0.01Hz
        @test abs(command_errors[switch+1]) < 0.01
        @test maximum(abs, command_errors[switch:end]) < 0.05
    end
end

@testset "Initial frequency and frequency rate" begin
    f0, a, T, N = 100.0, 50.0, 1e-3, 1000
    Δt = T * s
    # The default states are zero and keep their types when propagated.
    for lf in (FirstOrderLF(), SecondOrderBilinearLF(), SecondOrderBoxcarLF(),
        ThirdOrderBilinearLF(), ThirdOrderBoxcarLF(), ThirdOrderAssistedBilinearLF())
        @test predict_frequency(lf, 0.0s) === 0.0Hz
        input = is_assisted(lf) ? (0.1, 1.0Hz) : 0.1
        bandwidth = is_assisted(lf) ? (10Hz, 2Hz) : 10Hz
        @test typeof(propagate(lf, input, Δt, bandwidth)) == typeof(lf)
        @test typeof(@inferred reinit(lf)) == typeof(lf)
        @test typeof(@inferred reinit(typeof(lf))) == typeof(lf)
    end

    # A third order filter started at the Doppler f0 + a·t of the signal, with the NCO at
    # its prediction for the first integration, sees no phase error and commands the
    # Doppler at the centre of each next integration.
    for (F, bandwidth) in ((ThirdOrderBilinearLF, 10Hz), (ThirdOrderBoxcarLF, 10Hz),
        (ThirdOrderAssistedBilinearLF, (10Hz, 2Hz)))
        lf = @inferred reinit(F; frequency = f0 * Hz, frequency_rate = a * Hz / s)
        @test lf == F(; frequency = f0 * Hz, frequency_rate = a * Hz / s)
        @test lf == F(f0 * Hz, a * Hz^2)
        @test get_frequency_rate(lf) ≈ a * Hz / s
        for τ in (0.0, T / 2, 0.1)
            @test predict_frequency(lf, τ * s) ≈ (f0 + a * τ) * Hz
        end
        _, δθs, commands = track_ramp(lf, bandwidth; f0, a, T, N,
            command = Float64(predict_frequency(lf, Δt / 2) / Hz), φ0 = a * T^2 / 12)
        @test maximum(abs, δθs) < 1e-12
        @test commands ≈ [f0 + a * (k + 1.5) * T for k in 0:N-1]
    end

    # Below third order: a constant Doppler. reinit silently ignores a frequency rate,
    # and the constructors take none.
    for (F, bandwidth) in ((FirstOrderLF, 2Hz), (SecondOrderBilinearLF, 5Hz), (SecondOrderBoxcarLF, 5Hz))
        lf = @inferred F(; frequency = f0 * Hz)
        @test lf == @test_logs reinit(F; frequency = f0 * Hz, frequency_rate = a * Hz / s)
        @test_throws MethodError F(; frequency = f0 * Hz, frequency_rate = a * Hz / s)
        @test predict_frequency(lf, 0.1s) == f0 * Hz
        _, δθs, commands = track_ramp(lf, bandwidth; f0, a = 0.0, T, N, command = f0)
        @test all(iszero, δθs)
        @test all(==(f0), commands)
    end

    # Units are converted into those of the filter states.
    lf = ThirdOrderBilinearLF(; frequency = 1u"kHz", frequency_rate = 1u"kHz/s")
    @test lf isa typeof(ThirdOrderBilinearLF())
    @test lf.x1 ≈ 1000.0Hz && lf.x2 ≈ 1000.0Hz^2

    # reinit keeps the element types and units of a filter.
    for (lf, has_rate) in (
        (ThirdOrderBilinearLF(0.0f0Hz, 0.0f0Hz^2), true),
        (ThirdOrderBoxcarLF(0.0Hz, 0.0Hz / s), true),
        (ThirdOrderAssistedBilinearLF(0.0f0Hz, 0.0f0Hz^2, 0.0f0u"ms"), true),
        (SecondOrderBoxcarLF(0.0f0u"kHz"), false),
        (FirstOrderLF(0.0f0Hz), false),
    )
        moved = @inferred reinit(lf; frequency = 1u"kHz", frequency_rate = has_rate ? 5Hz / s : 0Hz / s)
        @test typeof(moved) == typeof(lf)
        @test predict_frequency(moved, 0.0s) ≈ 1000Hz
        @test get_frequency_rate(moved) ≈ (has_rate ? 5Hz / s : 0Hz / s)
    end

    # Reinitialising from a filter's own estimates moves them, e.g. onto another
    # reference frequency.
    lf, = track_ramp(ThirdOrderBoxcarLF(), 10Hz; f0, a, T, N = 5000)
    moved = @inferred reinit(lf; frequency = predict_frequency(lf, 0.0s) - 30Hz,
        frequency_rate = get_frequency_rate(lf) - 5Hz / s)
    for τ in (-T / 2, 0.0, 0.1)
        @test predict_frequency(moved, τ * s) ≈
            predict_frequency(lf, τ * s) - 30Hz - 5Hz / s * τ * s
    end
    @test filter_loop(moved, 0.0, Δt, 10Hz)[1] ≈
        filter_loop(lf, 0.0, Δt, 10Hz)[1] - 30Hz - 5Hz / s * 3Δt / 2
end
