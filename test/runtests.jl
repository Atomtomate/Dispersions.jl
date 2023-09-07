using Test

using Dispersions
using Random


grid_list_D = [2, 3, 3, 2, 2, 2, 2, 2, 3]
grid_list = ["2Dsc-1.3", "3Dsc-1.3", "fcc-1.4", "p6m-1.5", "p6m-1.5-1.7-0.0", "2Dsc--1.6", "2Dsc-1.7--1.8-1.9", "2Dsc-1.7--1.8--1.9", "bcc-1.1"]
num_eps = 1e-8
rng = MersenneTwister(0);

include("./helper_functions.jl")
@testset "cP" begin
    include("./Types.jl")
    include("./cP.jl")
    #include("./IO_SC.jl")    
end

@testset "KGrid" begin
    include("./KGrid.jl")
end

@testset "cPnn" begin
    include("./cPnn.jl")
end

@testset "cF" begin
    include("./cF.jl")
end

@testset "cI" begin
    include("./cI.jl")
end

@testset "hexagonal" begin
    include("./hexagonal.jl")
end


@testset "BZ Integration" begin
    include("./BZIntegration.jl")
end

@testset "common" begin
    include("./common.jl")
end
