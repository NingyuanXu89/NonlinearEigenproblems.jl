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

function region_corners(region)
    half_width = isa(region, SquareRegion) ? region.half_side : region.half_width
    half_height = isa(region, SquareRegion) ? region.half_side : region.half_height
    return [region.center + sx*half_width + im*sy*half_height
            for sx in (-1, 1), sy in (-1, 1)]
end

function inside_contour(contour, z)
    dz = z - contour.center
    return (real(dz) / contour.radius[1])^2 + (imag(dz) / contour.radius[2])^2 <=
        1 + 10eps()
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
    @test ellipse.radius == (sqrt(2) * rect.half_width, sqrt(2) * rect.half_height)
    @test all(corner -> inside_contour(ellipse, corner), region_corners(rect))
    square_ellipse = enclosing_contour(square; shape=:ellipse)
    @test square_ellipse.radius == circle.radius
    @test all(corner -> inside_contour(square_ellipse, corner), region_corners(square))
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

    square_boundary = region_boundary(square; n=2)
    @test square_boundary.z == [0.5 + 1.5im, 1.5 + 1.5im, 1.5 + 2.5im,
                               0.5 + 2.5im, 0.5 + 1.5im]
    @test square_boundary.x == real.(square_boundary.z)
    @test square_boundary.y == imag.(square_boundary.z)
    @test first(square_boundary.z) == last(square_boundary.z)
    @test_throws ErrorException region_boundary(square; n=1)

    rect_boundary = region_boundary(rect; n=3)
    @test length(rect_boundary.z) == 9
    @test extrema(rect_boundary.x) == (-1.0, 3.0)
    @test extrema(rect_boundary.y) == (-1.25, -0.75)

    contour = EllipseContour(1 + 2im, (2.0, 1.0))
    contour_closed = contour_boundary(contour; n=4)
    @test length(contour_closed.z) == 5
    @test contour_closed.z ≈ [3 + 2im, 1 + 3im, -1 + 2im, 1 + 1im, 3 + 2im]
    @test first(contour_closed.z) ≈ last(contour_closed.z)

    contour_open = contour_boundary(contour; n=4, closed=false)
    @test length(contour_open.z) == 4
    @test first(contour_open.z) != last(contour_open.z)
    @test contour_boundary(1 + 2im, (2.0, 1.0); n=4).z ≈ contour_closed.z
    @test contour_boundary(; σ=1 + 2im, radius=(2.0, 1.0), n=4).z ≈ contour_closed.z
    @test_throws ErrorException contour_boundary(contour; n=0)

    depth_records = collect_region_boundaries(square; depth=1)
    @test length(depth_records) == 5
    @test getproperty.(depth_records, :id) == ["root", "1", "2", "3", "4"]
    @test getproperty.(depth_records, :level) == [0, 1, 1, 1, 1]
    @test depth_records[1].parent_id === nothing
    @test depth_records[2].parent_id == "root"
    @test depth_records[2].child_index == 1
    @test depth_records[2].region.center == 1.25 + 2.25im

    selected_records = collect_region_boundaries(square; selected=(1, 4))
    @test getproperty.(selected_records, :id) == ["root", "1", "1.4"]
    @test selected_records[3].parent_id == "1"
    @test selected_records[3].child_index == 4
    @test selected_records[3].region.center == 1.125 + 2.125im
    @test getproperty.(collect_region_boundaries(square; selected=[(1,), (2, 3)]), :id) ==
        ["root", "1", "2", "2.3"]
    @test_throws ErrorException collect_region_boundaries(square; depth=-1)
    @test_throws ErrorException collect_region_boundaries(square; selected=(5,))

    default_contours = collect_contour_boundaries(rect; n=4)
    @test length(default_contours) == 1
    @test default_contours[1].id == "root"
    @test default_contours[1].contour.radius[1] ≈ sqrt(rect.half_width^2 + rect.half_height^2)
    @test default_contours[1].contour.radius[1] == default_contours[1].contour.radius[2]
    @test length(default_contours[1].z) == 5

    ellipse_contours = collect_contour_boundaries(rect; n=4, shape=:ellipse)
    @test ellipse_contours[1].contour.radius ==
        (sqrt(2) * rect.half_width, sqrt(2) * rect.half_height)
    @test all(corner -> inside_contour(ellipse_contours[1].contour, corner),
              region_corners(rect))
end

@testset "SIM indicator and screening" begin
    λ0 = 0.2 + 0.1im
    nep = diagonal_linear_pep([λ0])

    active_region = SquareRegion(λ0, 0.1)
    inactive_region = SquareRegion(2 + 0im, 0.1)

    @test sim_indicator(nep, active_region, N=16) ≈ 1 atol=1e-12
    @test sim_indicator(nep, active_region, N=16, method=:xisun) ==
        sim_indicator(nep, active_region, N=16)
    @test sim_indicator(nep, inactive_region, N=16) < 0.1

    active = sim_screen(nep, active_region, N=16, threshold=0.1)
    active_xisun = sim_screen(nep, active_region, N=16, threshold=0.1,
                              method=:xisun)
    inactive = sim_screen(nep, inactive_region, N=16, threshold=0.1)
    @test active.region == active_region
    @test active.active
    @test inactive.region == inactive_region
    @test !inactive.active
    @test active.threshold == 0.1
    @test active.N == 16
    @test active.indicator == sim_indicator(nep, active_region, N=16)
    @test active_xisun.indicator == active.indicator
    @test active_xisun.active == active.active
    @test active_xisun.threshold == active.threshold
    @test active_xisun.N == active.N
    @test active_xisun.probe_norm == active.probe_norm
    @test_throws MethodError sim_indicator(nep, active_region, N=16, threshold=0.1)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:unknown)
    @test_throws ErrorException sim_screen(nep, active_region, N=16,
                                           method=:unknown)

    rect = RectangularRegion(λ0, 0.2, 0.05)
    @test sim_screen(nep, rect, N=16).active

    results = sim_screen_regions(nep, [active_region, inactive_region], N=16)
    @test map(r -> r.active, results) == [true, false]
end

@testset "SIM diagnostics plumbing" begin
    nep = diagonal_linear_pep([0.2 + 0.1im])
    region = SquareRegion(0.2 + 0.1im, 0.1)
    contour = enclosing_contour(region; shape=:circle)

    result = sim_screen(nep, region, N=16, seed=99, method=:xisun)
    @test result.diagnostics.method == :xisun
    @test result.diagnostics.K == 0
    @test result.diagnostics.k == 1
    @test isempty(result.diagnostics.moment_strengths)
    @test isempty(result.diagnostics.moment_norms)
    @test isempty(result.diagnostics.moment_singular_values)
    @test result.diagnostics.seed == 99
    @test result.diagnostics.contour_center == contour.center
    @test result.diagnostics.contour_radius == contour.radius

    legacy = SIMResult(region, result.indicator, result.active, result.threshold,
                       result.N, result.probe_norm)
    @test legacy.indicator == result.indicator
    @test legacy.active == result.active
    @test legacy.diagnostics.method == :xisun
    @test legacy.diagnostics.seed == 10
end

@testset "SIM moment_norm screening" begin
    λ0 = 0.2 + 0.1im
    nep = diagonal_linear_pep([λ0, -0.3 + 0.2im])
    active_region = SquareRegion(λ0, 0.1)
    inactive_region = SquareRegion(2 + 0im, 0.1)
    contour = enclosing_contour(active_region; shape=:circle)

    result = sim_screen(nep, active_region, N=16, method=:moment_norm,
                        k=2, K=2, seed=11, threshold=0.1)
    @test result isa SIMResult
    @test result.diagnostics.method == :moment_norm
    @test result.diagnostics.K == 2
    @test result.diagnostics.k == 2
    @test length(result.diagnostics.moment_strengths) == 3
    @test length(result.diagnostics.moment_norms) == 3
    @test length(result.diagnostics.moment_singular_values) == 3
    @test all(isfinite, result.diagnostics.moment_strengths)
    @test all(isfinite, result.diagnostics.moment_norms)
    @test all(result.diagnostics.moment_strengths .>= 0)
    @test all(result.diagnostics.moment_norms .>= 0)
    @test result.diagnostics.seed == 11
    @test result.diagnostics.contour_center == contour.center
    @test result.diagnostics.contour_radius == contour.radius
    @test result.indicator == maximum(result.diagnostics.moment_strengths)

    @test sim_indicator(nep, active_region, N=16, method=:moment_norm,
                        k=2, K=2, seed=11) == result.indicator

    repeated = sim_screen(nep, active_region, N=16, method=:moment_norm,
                          k=2, K=2, seed=11, threshold=0.1)
    changed = sim_screen(nep, active_region, N=16, method=:moment_norm,
                         k=2, K=2, seed=12, threshold=0.1)
    @test repeated.diagnostics.moment_strengths == result.diagnostics.moment_strengths
    @test repeated.diagnostics.moment_norms == result.diagnostics.moment_norms
    @test repeated.diagnostics.moment_singular_values ==
        result.diagnostics.moment_singular_values
    @test changed.diagnostics.moment_strengths != result.diagnostics.moment_strengths

    Random.seed!(1234)
    expected_first = rand()
    expected_second = rand()
    Random.seed!(1234)
    @test rand() == expected_first
    sim_screen(nep, active_region, N=16, method=:moment_norm,
               k=2, K=2, seed=99)
    @test rand() == expected_second

    active_indicator = sim_indicator(nep, active_region, N=16,
                                     method=:moment_norm, k=2, K=2, seed=21)
    inactive_indicator = sim_indicator(nep, inactive_region, N=16,
                                       method=:moment_norm, k=2, K=2, seed=21)
    @test active_indicator > inactive_indicator
    threshold = (active_indicator + inactive_indicator) / 2
    @test sim_screen(nep, active_region, N=16, method=:moment_norm,
                     k=2, K=2, seed=21, threshold=threshold).active
    @test !sim_screen(nep, inactive_region, N=16, method=:moment_norm,
                      k=2, K=2, seed=21, threshold=threshold).active

    explicit_vector = ComplexF64[2, 0]
    vector_result = sim_screen(nep, active_region, N=16, method=:moment_norm,
                               k=1, K=1, probe=explicit_vector)
    @test vector_result.probe_norm == norm(explicit_vector)

    explicit_matrix = ComplexF64[1 0; 0 1]
    matrix_result = sim_screen(nep, active_region, N=16, method=:moment_norm,
                               k=2, K=1, probe=explicit_matrix)
    @test matrix_result.probe_norm == norm(explicit_matrix)

    @test_throws ErrorException sim_indicator(nep, active_region, N=0,
                                              method=:moment_norm)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:moment_norm, k=0)
    @test_throws ErrorException sim_indicator(nep, active_region, N=1,
                                              method=:moment_norm, k=2)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:moment_norm, K=-1)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:moment_norm, K=4)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:moment_norm, k=2,
                                              probe=explicit_vector)
    @test_throws ErrorException sim_indicator(nep, active_region, N=16,
                                              method=:moment_norm, k=2,
                                              probe=ones(ComplexF64, 3, 2))
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
    @test skipped.reason == :legacy_xisun_inactive
    @test skipped.recommended_kwargs == NamedTuple()

    no_info = sim_contour_decision(active)
    @test no_info.action == :run_beyn
    @test no_info.reason == :active

    nep = diagonal_linear_pep([0.2 + 0.1im, -0.3 + 0.2im])
    region = SquareRegion(0.2 + 0.1im, 0.1)
    moment_indicator = sim_indicator(nep, region, N=16, method=:moment_norm,
                                     k=2, K=2, seed=3)

    active_moment = sim_screen(nep, region, N=16, method=:moment_norm,
                               k=2, K=2, seed=3,
                               threshold=moment_indicator / 2)
    moment_no_info = sim_contour_decision(active_moment)
    @test moment_no_info.action == :run_beyn
    @test moment_no_info.reason == :moment_detected

    far_inactive_moment = sim_screen(nep, region, N=16, method=:moment_norm,
                                     k=2, K=2, seed=3,
                                     threshold=20 * moment_indicator)
    far_decision = sim_contour_decision(far_inactive_moment)
    @test far_decision.action == :skip
    @test far_decision.reason == :moment_norm_inactive

    borderline_moment = sim_screen(nep, region, N=16, method=:moment_norm,
                                   k=2, K=2, seed=3,
                                   threshold=5 * moment_indicator)
    borderline_decision = sim_contour_decision(borderline_moment)
    @test borderline_decision.action == :subdivide_or_rerun
    @test borderline_decision.reason == :weak_moment_strength

    no_longer_borderline = sim_contour_decision(borderline_moment;
                                                borderline_factor=2)
    @test no_longer_borderline.action == :skip
    @test no_longer_borderline.reason == :moment_norm_inactive

    moment_within_target_info = (; estimated_rank=3, capacity=10, number_returned=3,
                                 residuals=[1e-12], singular_values=[1.0, 1e-2, 1e-12],
                                 rank_drop_tol=1e-8)
    within_target = sim_contour_decision(active_moment, moment_within_target_info;
                                         target_capacity=3, residual_tol=1e-8)
    @test within_target.action == :accept_region
    @test within_target.reason == :moment_rank_within_target

    above_target_info = (; estimated_rank=4, capacity=10, number_returned=4,
                         residuals=[1e-12], singular_values=[1.0, 1e-2, 1e-3, 1e-12],
                         rank_drop_tol=1e-8)
    above_target = sim_contour_decision(active_moment, above_target_info;
                                        target_capacity=3, residual_tol=1e-8)
    @test above_target.action == :subdivide_or_rerun
    @test above_target.reason == :moment_rank_above_target

    moment_no_pairs = (; estimated_rank=3, capacity=10, number_returned=0,
                       residuals=Float64[], singular_values=[1.0, 1e-2, 1e-12],
                       rank_drop_tol=1e-8)
    @test sim_contour_decision(active_moment, moment_no_pairs;
                               target_capacity=3).reason == :active_no_eigenpairs

    moment_poor_residuals = (; estimated_rank=3, capacity=10, number_returned=3,
                             residuals=[1e-2],
                             singular_values=[1.0, 1e-2, 1e-12],
                             rank_drop_tol=1e-8)
    @test sim_contour_decision(active_moment, moment_poor_residuals;
                               target_capacity=3, residual_tol=1e-8).reason ==
        :poor_residuals

    moment_unclear_rank = (; estimated_rank=1, capacity=10, number_returned=1,
                           residuals=[1e-12], singular_values=[1.0, 1e-8],
                           rank_drop_tol=1e-8)
    @test sim_contour_decision(active_moment, moment_unclear_rank;
                               target_capacity=3).reason == :unclear_rank_gap

    @test_throws ErrorException sim_contour_decision(active_moment,
                                                     moment_within_target_info;
                                                     target_capacity=0)

    low_rank = (; estimated_rank=1, capacity=4, number_returned=1,
                residuals=[1e-12], singular_values=[1.0, 1e-12, 1e-14],
                rank_drop_tol=1e-8)
    accepted = sim_contour_decision(active, low_rank; residual_tol=1e-8)
    @test accepted.action == :accept_beyn
    @test accepted.reason == :low_estimated_rank

    xisun_with_target = sim_contour_decision(active, low_rank;
                                             residual_tol=1e-8,
                                             target_capacity=1)
    @test xisun_with_target.action == accepted.action
    @test xisun_with_target.reason == accepted.reason

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

@testset "SIM auto-subdivision helper" begin
    λs = [-0.5 - 0.5im, -0.5 + 0.5im, 0.5 - 0.5im, 0.5 + 0.5im]
    nep = diagonal_linear_pep(λs)
    root = SquareRegion(0, 1.0)
    empty_region = SquareRegion(3 + 0im, 0.2)
    screen_N = 16
    screen_k = 1
    screen_K = 1
    target_capacity = screen_k * screen_K
    verify_kwargs = (; k=4, K=1, N=64, rank_drop_tol=1e-8)

    regions = sim_subdivide_active_regions(nep, root;
                                           N=screen_N, k=screen_k, K=screen_K,
                                           threshold=1e-8, seed=4,
                                           max_depth=2,
                                           target_capacity=target_capacity,
                                           verify_kwargs=verify_kwargs)
    @test !isempty(regions)
    @test all(region -> region isa Union{SquareRegion,RectangularRegion}, regions)
    @test all(λ -> any(region -> inside_region(region, λ), regions), λs)

    traced = sim_subdivide_active_regions(nep, root;
                                          N=screen_N, k=screen_k, K=screen_K,
                                          threshold=1e-8, seed=4,
                                          max_depth=2,
                                          target_capacity=target_capacity,
                                          verify_kwargs=verify_kwargs,
                                          return_trace=true)
    @test traced.regions == regions
    @test !isempty(traced.trace)
    @test all(entry -> hasproperty(entry, :path), traced.trace)
    @test all(entry -> hasproperty(entry, :status), traced.trace)
    @test all(entry -> hasproperty(entry, :first_decision), traced.trace)
    @test all(entry -> hasproperty(entry, :final_decision), traced.trace)
    accepted_trace = filter(entry -> entry.status == :accepted, traced.trace)
    accepted_paths = getproperty.(accepted_trace, :path)
    accepted_boundaries = collect_region_boundaries(root; selected=accepted_paths)
    @test getproperty.(accepted_trace, :region) == traced.regions
    @test all(path -> path in getproperty.(accepted_boundaries, :path), accepted_paths)
    @test :subdivided_over_target in getproperty.(traced.trace, :status)
    @test :accepted in getproperty.(traced.trace, :status)

    for region in regions
        contour = enclosing_contour(region; shape=:circle)
        info = contour_block_SS_info(nep; contour_parameters(contour)...,
                                     seed=4, verify_kwargs...)
        sim = sim_screen(nep, region; method=:moment_norm,
                         N=screen_N, k=screen_k, K=screen_K,
                         threshold=1e-8, seed=4)
        decision = sim_contour_decision(sim, info;
                                        target_capacity=target_capacity)
        @test decision.action == :accept_region
        @test info.estimated_rank <= target_capacity
    end

    skipped = sim_subdivide_active_regions(nep, empty_region;
                                           N=screen_N, k=screen_k, K=screen_K,
                                           threshold=1e-8, seed=4,
                                           max_depth=2,
                                           target_capacity=target_capacity,
                                           verify_kwargs=verify_kwargs)
    @test isempty(skipped)

    skipped_trace = sim_subdivide_active_regions(nep, empty_region;
                                                 N=screen_N, k=screen_k,
                                                 K=screen_K,
                                                 threshold=1e-8, seed=4,
                                                 max_depth=2,
                                                 target_capacity=target_capacity,
                                                 verify_kwargs=verify_kwargs,
                                                 return_trace=true)
    @test isempty(skipped_trace.regions)
    @test getproperty.(skipped_trace.trace, :status) == [:skipped_inactive]

    @test_throws ErrorException sim_subdivide_active_regions(nep, root;
                                                            N=screen_N,
                                                            k=screen_k,
                                                            K=screen_K,
                                                            threshold=1e-8,
                                                            seed=4,
                                                            max_depth=0,
                                                            target_capacity=target_capacity,
                                                            verify_kwargs=verify_kwargs)

    @test_logs (:warn, r"omitting unresolved SIM subdivision region") begin
        omitted = sim_subdivide_active_regions(nep, root;
                                               N=screen_N, k=screen_k,
                                               K=screen_K,
                                               threshold=1e-8, seed=4,
                                               max_depth=0,
                                               target_capacity=target_capacity,
                                               verify_kwargs=verify_kwargs,
                                               strict=false)
        @test isempty(omitted)
    end

    @test_logs (:warn, r"omitting unresolved SIM subdivision region") begin
        omitted_trace = sim_subdivide_active_regions(nep, root;
                                                     N=screen_N, k=screen_k,
                                                     K=screen_K,
                                                     threshold=1e-8, seed=4,
                                                     max_depth=0,
                                                     target_capacity=target_capacity,
                                                     verify_kwargs=verify_kwargs,
                                                     strict=false,
                                                     return_trace=true)
        @test isempty(omitted_trace.regions)
        @test getproperty.(omitted_trace.trace, :status) == [:omitted_unresolved]
    end
end
