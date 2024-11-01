using Test

using Dispersions
using Random
using LinearAlgebra


grid_list = ["2Dsc-1.3", "3Dsc-1.3", "fcc-1.4", "2Dsc--1.6", "2Dsc-1.7--1.8-1.9", "2Dsc-1.7--1.8--1.9", "bcc-1.1"]
num_eps = 1e-8
rng = MersenneTwister(0);

include("helper_functions.jl")
@testset "cP" begin
    include("Types.jl")
    include("cP.jl")
    #include("./IO_SC.jl")    
end

@testset "KGrid" begin
    include("KGrid.jl")
end

@testset "cPnn" begin
    include("cPnn.jl")
end

@testset "cF" begin
    include("cF.jl")
end

@testset "cI" begin
    include("cI.jl")
end


@testset "BZ Integration" begin
    include("BZIntegration.jl")
end

@testset "common" begin
    include("common.jl")
end

@testset "helpers" begin
    include("helpers.jl")
end
