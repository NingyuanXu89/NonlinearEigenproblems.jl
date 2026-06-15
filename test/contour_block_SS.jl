# Run tests on the block SS contour integral method

using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random


@testset "block SS" begin
    nep=nep_gallery("dep0",3)

    # Circle
    λ,V=contour_block_SS(nep,radius=1.0,N=1000,σ=0.1,k=3,K=3);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

    # Ellipse
    λ,V=contour_block_SS(nep,radius=[1.0,2.0],N=1000,σ=0.1,k=3,K=3);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

    # JSIAM Mode
    λ,V=contour_block_SS(nep,radius=1.0,N=1000,σ=0.1,k=3,K=4,Shat_mode=:JSIAM);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

end

@testset "block SS diagnostics" begin
    nep = nep_gallery("dep0", 3)

    info = contour_block_SS_info(nep, radius=1.0, N=1000, σ=0.1, k=3, K=3)
    λ, V = contour_block_SS(nep, radius=1.0, N=1000, σ=0.1, k=3, K=3)

    @test info.lambda == λ
    @test info.V == V
    @test info.capacity == 9
    @test length(info.singular_values) == info.capacity
    @test info.estimated_rank <= info.capacity
    @test info.number_returned == length(info.lambda)
    @test length(info.residuals) == info.number_returned
    @test length(info.inside_contour) == info.number_returned
    @test info.estimated_rank ==
        count(info.singular_values / info.singular_values[1] .> info.rank_drop_tol)

    seeded_a = contour_block_SS_info(nep, radius=1.0, N=200, σ=0.1,
                                     k=2, K=2, seed=77)
    seeded_b = contour_block_SS_info(nep, radius=1.0, N=200, σ=0.1,
                                     k=2, K=2, seed=77)
    seeded_c = contour_block_SS_info(nep, radius=1.0, N=200, σ=0.1,
                                     k=2, K=2, seed=78)
    @test seeded_a.singular_values == seeded_b.singular_values
    @test seeded_a.singular_values != seeded_c.singular_values

    Random.seed!(1234)
    expected_first = rand()
    expected_second = rand()
    Random.seed!(1234)
    @test rand() == expected_first
    contour_block_SS_info(nep, radius=1.0, N=200, σ=0.1,
                          k=2, K=2, seed=99)
    @test rand() == expected_second
end
