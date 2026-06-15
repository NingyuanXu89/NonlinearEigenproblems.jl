using LinearAlgebra
using NonlinearEigenproblems

function diagonal_linear_nep(λ)
    n = length(λ)
    A0 = Matrix(Diagonal(-complex.(λ)))
    A1 = Matrix{ComplexF64}(I, n, n)
    return PEP([A0, A1])
end

function deduplicate_eigenpairs(λ, V; tol)
    keep = Int[]
    for i in eachindex(λ)
        if all(abs(λ[i] - λ[j]) > tol for j in keep)
            push!(keep, i)
        end
    end
    return λ[keep], V[:, keep]
end

function run_sim_workflow_demo()
    known_λ = ComplexF64[-0.35 + 0.10im, 0.20 + 0.05im, 0.70 - 0.20im]
    nep = diagonal_linear_nep(known_λ)

    region = SquareRegion(0.0 + 0.0im, 0.75)

    screen_N = 32
    screen_k = 1
    screen_K = 1
    threshold = 1e-8
    seed = 1
    max_depth = 3
    target_capacity = screen_k * screen_K

    additional_moments = 2
    recovery_K = screen_K + additional_moments
    recovery_N = 96
    contour_shape = :circle
    dedup_tol = 1e-8

    subdivision = sim_subdivide_active_regions(
        nep, region;
        N=screen_N,
        k=screen_k,
        K=screen_K,
        threshold=threshold,
        seed=seed,
        max_depth=max_depth,
        target_capacity=target_capacity,
        verify_kwargs=(; k=screen_k, K=recovery_K, N=recovery_N,
                       rank_drop_tol=1e-9),
        return_trace=true
    )

    regions = subdivision.regions
    trace = subdivision.trace

    @assert !isempty(regions)

    accepted_paths = getproperty.(filter(t -> t.status == :accepted, trace), :path)
    accepted_boundaries = collect_region_boundaries(region; selected=accepted_paths)
    tree_boundaries = collect_region_boundaries(region;
                                                selected=getproperty.(trace, :path))

    candidate_λ = ComplexF64[]
    candidate_V = Matrix{ComplexF64}(undef, size(nep, 1), 0)
    candidate_residuals = Float64[]
    recovery_ranks = Int[]

    for active_region in regions
        contour = enclosing_contour(active_region; shape=contour_shape)
        params = contour_parameters(contour)
        info = contour_block_SS_info(
            nep;
            params...,
            k=screen_k,
            K=recovery_K,
            N=recovery_N,
            seed=seed,
            rank_drop_tol=1e-9
        )
        push!(recovery_ranks, info.estimated_rank)

        inside = filter(i -> inside_region(active_region, info.lambda[i]),
                        eachindex(info.lambda))
        append!(candidate_λ, info.lambda[inside])
        candidate_V = hcat(candidate_V, info.V[:, inside])
        append!(candidate_residuals, info.residuals[inside])
    end

    refined_λ = ComplexF64[]
    refined_V = Matrix{ComplexF64}(undef, size(nep, 1), 0)
    for i in eachindex(candidate_λ)
        λ_refined, v_refined = newton(nep; λ=candidate_λ[i],
                                      v=candidate_V[:, i],
                                      tol=1e-12,
                                      maxit=5)
        push!(refined_λ, λ_refined)
        refined_V = hcat(refined_V, v_refined)
    end

    solved_λ, solved_V = deduplicate_eigenpairs(refined_λ, refined_V;
                                                tol=dedup_tol)

    sort_by_known_distance(λ) = sort(λ, by=z -> minimum(abs.(known_λ .- z)))
    sorted_solved_λ = sort_by_known_distance(solved_λ)

    @assert length(solved_λ) == 3
    @assert size(solved_V, 2) == 3
    @assert maximum(map(z -> minimum(abs.(known_λ .- z)), solved_λ)) < 1e-10
    @assert all(rank -> rank <= target_capacity, recovery_ranks)

    return (known_λ=known_λ,
            nep=nep,
            region=region,
            regions=regions,
            trace=trace,
            accepted_paths=accepted_paths,
            accepted_boundaries=accepted_boundaries,
            tree_boundaries=tree_boundaries,
            candidate_λ=candidate_λ,
            candidate_V=candidate_V,
            candidate_residuals=candidate_residuals,
            recovery_ranks=recovery_ranks,
            solved_λ=solved_λ,
            solved_V=solved_V,
            sorted_solved_λ=sorted_solved_λ)
end

demo = run_sim_workflow_demo()
known_λ = demo.known_λ
nep = demo.nep
region = demo.region
regions = demo.regions
trace = demo.trace
accepted_paths = demo.accepted_paths
accepted_boundaries = demo.accepted_boundaries
tree_boundaries = demo.tree_boundaries
candidate_λ = demo.candidate_λ
candidate_V = demo.candidate_V
candidate_residuals = demo.candidate_residuals
recovery_ranks = demo.recovery_ranks
solved_λ = demo.solved_λ
solved_V = demo.solved_V
sorted_solved_λ = demo.sorted_solved_λ

println("SIM moment-norm workflow demo")
println("accepted regions: ", length(regions))
println("trace entries: ", length(trace))
println("accepted region paths: ", accepted_paths)
println("candidate eigenvalues before Newton/de-dup: ", candidate_λ)
println("max contour residual: ", maximum(candidate_residuals))
println("solved eigenvalues: ", sorted_solved_λ)
