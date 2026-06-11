using Test

@testset "SIM workflow demo" begin
    include(joinpath(@__DIR__, "..", "docs", "src", "sim_workflow_demo.jl"))
    @test true
end
