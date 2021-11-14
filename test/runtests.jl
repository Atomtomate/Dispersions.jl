using Test

using Dispersions


grid_list = ["2Dsc-1.3", "3Dsc-1.3", "p6m-1.1"]

include("./helper_functions.jl")
@testset "SC" begin
    include("./Types.jl")
    include("./SC.jl")    
    #include("./IO_SC.jl")    
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
