using TargeneCore
using Test

TESTDIR = joinpath(pkgdir(TargeneCore), "test")
@testset "Test TargeneCore" begin
    @test include(joinpath(TESTDIR, "utils.jl"))
    @test include(joinpath(TESTDIR, "dataset.jl"))
    @test include(joinpath(TESTDIR, "plots.jl"))
    @test include(joinpath(TESTDIR, "confounders.jl"))
    @test include(joinpath(TESTDIR, "inputs_from_estimands.jl"))
    @test include(joinpath(TESTDIR, "inputs_from_config.jl"))
    @test include(joinpath(TESTDIR, "inputs_from_gwas_config.jl"))
end
