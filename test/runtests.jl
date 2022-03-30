using Test

using Dispersions
using Random


grid_list = ["2Dsc-1.3", "3Dsc-1.3", "fcc-1.4", "p6m-1.5"]
num_eps = 1e-8
rng = MersenneTwister(0);

include("./helper_functions.jl")
@testset "SC" begin
    include("./Types.jl")
    include("./SC.jl")    
    #include("./IO_SC.jl")    
end

@testset "FCC" begin
    include("./FCC.jl")
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
