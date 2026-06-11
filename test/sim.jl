using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random

function diagonal_linear_pep(λv)
    n = length(λv)
    A0 = Matrix(Diagonal(-complex.(λv)))
    A1 = Matrix{ComplexF64}(I, n, n)
    return PEP([A0, A1])
end

@testset "SIM region helpers" begin
    @test_throws ErrorException SquareRegion(0, 0)
    @test_throws ErrorException RectangularRegion(0, 1, 0)
    @test_throws ErrorException EllipseContour(0, (1, 0))

    square = SquareRegion(1 + 2im, 0.5)
    rect = RectangularRegion(1 - im, 2.0, 0.25)

    @test inside_region(square, 1.4 + 2.5im)
    @test !inside_region(square, 1.6 + 2im)
    @test inside_region(rect, 3 - 0.75im)
    @test !inside_region(rect, 3.1 - im)

    circle = enclosing_contour(square)
    @test circle.center == square.center
    @test circle.radius[1] ≈ sqrt(2) * square.half_side
    @test circle.radius[1] == circle.radius[2]

    ellipse = enclosing_contour(rect; shape=:ellipse)
    @test ellipse.radius == (rect.half_width, rect.half_height)
    params = contour_parameters(ellipse)
    @test params.σ == rect.center
    @test params.radius == ellipse.radius

    square_children = subdivide(square)
    @test length(square_children) == 4
    @test all(child -> child.half_side == square.half_side / 2, square_children)
    @test all(child -> inside_region(square, child.center), square_children)

    rect_children = subdivide(rect)
    @test length(rect_children) == 4
    @test all(child -> child.half_width == rect.half_width / 2, rect_children)
    @test all(child -> child.half_height == rect.half_height / 2, rect_children)
    @test all(child -> inside_region(rect, child.center), rect_children)
end

@testset "SIM indicator and screening" begin
    λ0 = 0.2 + 0.1im
    nep = diagonal_linear_pep([λ0])

    active_region = SquareRegion(λ0, 0.1)
    inactive_region = SquareRegion(2 + 0im, 0.1)

    @test sim_indicator(nep, active_region, N=16) ≈ 1 atol=1e-12
    @test sim_indicator(nep, inactive_region, N=16) < 0.1

    active = sim_screen(nep, active_region, N=16, threshold=0.1)
    inactive = sim_screen(nep, inactive_region, N=16, threshold=0.1)
    @test active.region == active_region
    @test active.active
    @test inactive.region == inactive_region
    @test !inactive.active
    @test active.threshold == 0.1
    @test active.N == 16
    @test active.indicator == sim_indicator(nep, active_region, N=16)
    @test_throws MethodError sim_indicator(nep, active_region, N=16, threshold=0.1)

    rect = RectangularRegion(λ0, 0.2, 0.05)
    @test sim_screen(nep, rect, N=16).active

    results = sim_screen_regions(nep, [active_region, inactive_region], N=16)
    @test map(r -> r.active, results) == [true, false]
end

@testset "SIM probe and validation" begin
    nep = diagonal_linear_pep([0.2 + 0.1im, -0.3 + 0.2im])
    region = SquareRegion(0.2 + 0.1im, 0.1)

    first = sim_screen(nep, region, N=16)
    second = sim_screen(nep, region, N=16)
    @test first.indicator == second.indicator
    @test first.probe_norm == second.probe_norm

    explicit_probe = ComplexF64[2, 0]
    explicit = sim_screen(nep, region, N=16, probe=explicit_probe)
    @test explicit.probe_norm == norm(explicit_probe)
    @test explicit.indicator ≈ sim_indicator(nep, region, N=16, probe=explicit_probe)

    Random.seed!(1234)
    expected_first = rand()
    expected_second = rand()
    Random.seed!(1234)
    @test rand() == expected_first
    sim_screen(nep, region, N=16, seed=99)
    @test rand() == expected_second

    @test_throws ErrorException sim_indicator(nep, region, N=15)
    @test_throws ErrorException sim_indicator(nep, region, N=0)
    @test_throws ErrorException sim_indicator(nep, region, N=16, probe=ComplexF64[1])
    @test_throws ErrorException sim_indicator(nep, region, N=16, probe=zeros(ComplexF64, 2))
end

@testset "SIM contour decision policy" begin
    inactive = SIMResult(SquareRegion(0, 1), 0.0, false, 0.1, 16, 1.0)
    active = SIMResult(SquareRegion(0, 1), 1.0, true, 0.1, 16, 1.0)

    skipped = sim_contour_decision(inactive)
    @test skipped.action == :skip
    @test skipped.reason == :inactive
    @test skipped.recommended_kwargs == NamedTuple()

    no_info = sim_contour_decision(active)
    @test no_info.action == :run_beyn
    @test no_info.reason == :active

    low_rank = (; estimated_rank=1, capacity=4, number_returned=1,
                residuals=[1e-12], singular_values=[1.0, 1e-12, 1e-14],
                rank_drop_tol=1e-8)
    accepted = sim_contour_decision(active, low_rank; residual_tol=1e-8)
    @test accepted.action == :accept_beyn
    @test accepted.reason == :low_estimated_rank

    near_capacity = (; estimated_rank=3, capacity=4, number_returned=3,
                     residuals=[1e-12], singular_values=[1.0, 1e-2, 1e-4, 1e-12],
                     rank_drop_tol=1e-8)
    rerouted = sim_contour_decision(active, near_capacity)
    @test rerouted.action == :run_block_SS_or_subdivide
    @test rerouted.reason == :near_capacity_rank

    no_pairs = (; estimated_rank=1, capacity=4, number_returned=0,
                residuals=Float64[], singular_values=[1.0, 1e-12],
                rank_drop_tol=1e-8)
    @test sim_contour_decision(active, no_pairs).reason == :active_no_eigenpairs

    poor_residuals = (; estimated_rank=1, capacity=4, number_returned=1,
                      residuals=[1e-2], singular_values=[1.0, 1e-12],
                      rank_drop_tol=1e-8)
    @test sim_contour_decision(active, poor_residuals; residual_tol=1e-8).reason ==
        :poor_residuals

    unclear_by_next_singular_value = (; estimated_rank=1, capacity=4, number_returned=1,
                                      residuals=[1e-12], singular_values=[1.0, 1e-8],
                                      rank_drop_tol=1e-8)
    unclear = sim_contour_decision(active, unclear_by_next_singular_value)
    @test unclear.action == :subdivide_or_rerun
    @test unclear.reason == :unclear_rank_gap

    unclear_by_ratio = (; estimated_rank=1, capacity=4, number_returned=1,
                        residuals=[1e-12], singular_values=[1.0, 0.2],
                        rank_drop_tol=10.0)
    @test sim_contour_decision(active, unclear_by_ratio).reason == :unclear_rank_gap

    missing_optional_fields = (; estimated_rank=1, capacity=4, number_returned=1)
    @test sim_contour_decision(active, missing_optional_fields).action == :accept_beyn

    missing_number_returned = (; estimated_rank=1, capacity=4,
                               residuals=[1e-12],
                               singular_values=[1.0, 1e-12],
                               rank_drop_tol=1e-8)
    missing_returned_decision = sim_contour_decision(active, missing_number_returned)
    @test missing_returned_decision.action == :accept_beyn
    @test missing_returned_decision.reason != :active_no_eigenpairs
end
