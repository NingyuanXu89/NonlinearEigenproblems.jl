# Run tests on Beyns contour integral method

using NonlinearEigenproblems
using Test
using Random
using LinearAlgebra


@testset "Beyn contour" begin

    nep=nep_gallery("dep0")
    @bench @testset "disk at origin" begin
        Random.seed!(0);

        λ,v=contour_beyn(nep,logger=displaylevel,radius=1,neigs=1,sanity_check=false)

        for i = 1:size(λ,1)
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*1000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*500
        end


        λ,v=contour_beyn(nep,logger=displaylevel,radius=0.4,k=2,N=100,sanity_check=false)
        M=compute_Mder(nep,λ[1])
        minimum(svdvals(M))
        @test minimum(svdvals(M))<eps()*1000

    end
    @bench @testset "shifted disk" begin

        λ,v=contour_beyn(nep,logger=displaylevel,σ=0.2,radius=1.0,
                         neigs=4,sanity_check=false)

        @test size(λ,1)==3
        for i = 1:3
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*10000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*10000
        end

    end

end

@testset "Beyn contour diagnostics" begin
    nep = nep_gallery("dep0")

    info = contour_beyn_info(nep, logger=displaylevel, radius=1.0,
                             neigs=1, k=2, N=100, sanity_check=false)
    λ, V = contour_beyn(nep, logger=displaylevel, radius=1.0,
                        neigs=1, k=2, N=100, sanity_check=false)

    @test info.lambda == λ
    @test info.V == V
    @test length(info.singular_values) == info.capacity
    @test info.estimated_rank <= info.capacity
    @test info.number_returned == length(info.lambda)
    @test length(info.residuals) == info.number_returned
    @test length(info.inside_contour) == info.number_returned
    @test info.estimated_rank ==
        count(info.singular_values / info.singular_values[1] .> info.rank_drop_tol)

    info_checked = contour_beyn_info(nep, logger=displaylevel, σ=0.2,
                                     radius=1.0, neigs=3, sanity_check=true)
    @test info_checked.number_returned == length(info_checked.lambda)
    @test length(info_checked.residuals) == info_checked.number_returned
    @test length(info_checked.inside_contour) == info_checked.number_returned
    @test all(info_checked.residuals .< sqrt(eps(Float64)))
end
