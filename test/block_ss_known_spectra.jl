using Test
using LinearAlgebra
using NonlinearEigenproblems

include(joinpath(@__DIR__, "nep_benchmarks", "block_ss_known_spectra_problems.jl"))

with_block_ss_kwargs(p) = merge(p, (; solver_kwargs=block_ss_solver_kwargs(p)))

clean_cases = with_block_ss_kwargs.([
    diagonal_polynomial_pep(),
    diagonal_trig_polynomial_nep(),
    shifted_diagonal_exponential_nep(),
    diagonal_delay_lambertw_nep(),
    similarity_diagonal_nep(),
])

stress_cases = with_block_ss_kwargs.([
    ill_scaled_diagonal_nep(),
    jordan_like_nep(),
    clustered_boundary_nep(),
])

@testset "block SS known spectra" begin
    for p in clean_cases
        @testset "$(p.name)" begin
            λ, V = contour_block_SS(p.nep; p.solver_kwargs...)

            @test all(inside_problem_contour.(Ref(p), λ))
            @test match_eigenvalues(λ, p.inside; atol=1e-8, rtol=1e-8).matched
            @test isempty(accepted_outside_values(λ, p.outside; atol=1e-8, rtol=1e-8))

            for i in eachindex(λ)
                @test norm(compute_Mlincomb(p.nep, λ[i], V[:, i])) / norm(V[:, i]) < 1e-8
            end
        end
    end
end

@testset "block SS known spectra diagnostics" begin
    for p in clean_cases
        @testset "$(p.name)" begin
            info = contour_block_SS_info(p.nep; p.solver_kwargs...)
            inside_λ = info.lambda[info.inside_contour]

            @test info.number_returned == length(info.lambda)
            @test length(info.singular_values) == info.capacity
            @test length(info.residuals) == info.number_returned
            @test length(info.inside_contour) == info.number_returned
            @test info.estimated_rank == length(p.inside)
            @test info.capacity == p.solver_kwargs.K * p.solver_kwargs.k
            @test all(info.residuals .< 1e-8)
            @test match_eigenvalues(inside_λ, p.inside; atol=1e-8, rtol=1e-8).matched
            @test isempty(accepted_outside_values(inside_λ, p.outside; atol=1e-8, rtol=1e-8))
        end
    end
end

@testset "block SS known spectra stress cases" begin
    for p in stress_cases
        @testset "$(p.name)" begin
            info = contour_block_SS_info(p.nep; p.solver_kwargs...)
            inside_λ = info.lambda[info.inside_contour]

            @test info.number_returned == length(info.lambda)
            @test length(info.singular_values) == info.capacity
            @test length(info.residuals) == info.number_returned
            @test length(info.inside_contour) == info.number_returned
            @test info.estimated_rank <= info.capacity
            @test info.capacity == p.solver_kwargs.K * p.solver_kwargs.k
            @test all(isfinite, info.lambda)
            @test all(isfinite, info.residuals)
            @test !isempty(inside_λ)
            @test all(inside_problem_contour.(Ref(p), inside_λ))

            if p.name == "ill_scaled_diagonal_nep"
                # This case intentionally varies diagonal component scales by orders of magnitude.
                @test all(info.residuals .< 1e-6)
                @test match_eigenvalues(inside_λ, p.inside; atol=1e-7, rtol=1e-7).matched
                @test isempty(accepted_outside_values(inside_λ, p.outside; atol=1e-7, rtol=1e-7))
            elseif p.name == "jordan_like_nep"
                @test all(info.residuals .< 1e-5)
                @test any(root -> match_eigenvalues(inside_λ, [root]; atol=1e-5, rtol=1e-5).matched,
                          p.inside)
                @test length(unique(p.inside)) < length(p.inside)
            elseif p.name == "clustered_boundary_nep"
                @test all(info.residuals .< 1e-5)
                @test any(root -> match_eigenvalues(inside_λ, [root]; atol=1e-5, rtol=1e-5).matched,
                          p.inside)
                @test any(root -> abs(abs(root - p.σ) - p.radius) < 0.05, p.inside)
            end
        end
    end
end
