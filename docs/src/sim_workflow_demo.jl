using LinearAlgebra
using NonlinearEigenproblems

function diagonal_linear_nep(λ)
    n = length(λ)
    A0 = Matrix(Diagonal(-complex.(λ)))
    A1 = Matrix{ComplexF64}(I, n, n)
    return PEP([A0, A1])
end

known_λ = ComplexF64[-0.35 + 0.10im, 0.20 + 0.05im, 0.70 - 0.20im]
nep = diagonal_linear_nep(known_λ)

region = SquareRegion(0.0 + 0.0im, 0.75)
regions = subdivide(region)
screened = sim_screen_regions(nep, regions; N=32, threshold=0.1, seed=1)
active = [result for result in screened if result.active]

@assert !isempty(active)

accepted_λ = ComplexF64[]
accepted_residuals = Float64[]
decisions = Symbol[]
refined_λ = ComplexF64[]

for result in active
    contour = enclosing_contour(result.region)
    params = contour_parameters(contour)
    info = contour_beyn_info(nep; params..., k=3, neigs=typemax(Int), N=64, tol=1e-10)
    decision = sim_contour_decision(result, info; residual_tol=1e-8)
    push!(decisions, decision.action)

    inside = filter(i -> inside_region(result.region, info.lambda[i]), eachindex(info.lambda))
    append!(accepted_λ, info.lambda[inside])
    append!(accepted_residuals, info.residuals[inside])

    # Optional postprocessing: refine accepted contour candidates with Newton.
    for i in inside
        λ_refined, _ = newton(nep; λ=info.lambda[i], v=info.V[:,i], tol=1e-12, maxit=5)
        push!(refined_λ, λ_refined)
    end
end

sort_by_known_distance(λ) = sort(λ, by=z -> minimum(abs.(known_λ .- z)))
accepted_λ = sort_by_known_distance(accepted_λ)
refined_λ = sort_by_known_distance(refined_λ)

@assert length(accepted_λ) == 3
@assert maximum(map(z -> minimum(abs.(known_λ .- z)), accepted_λ)) < 1e-8
@assert maximum(map(z -> minimum(abs.(known_λ .- z)), refined_λ)) < 1e-10

println("SIM workflow demo")
println("active regions: ", length(active), " / ", length(regions))
println("policy decisions: ", decisions)
println("accepted eigenvalues: ", accepted_λ)
println("max contour residual: ", maximum(accepted_residuals))
println("optional Newton-refined eigenvalues: ", refined_λ)
