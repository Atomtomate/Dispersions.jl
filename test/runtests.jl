using Test

using Dispersions


grid_list = ["2Dsc-1.3", "3Dsc-1.3","fcc-1.4"] # "p6m-1.1", 
grid_list_D = [2,3,3]

include("./helper_functions.jl")
@testset "IO" begin
    include("./IO.jl")
end

@testset "common" begin
    include("./common.jl")
end

@testset "KGrid" begin
    include("./Types.jl")
    include("./KGrid.jl")
end

@testset "BZ Integration" begin
    include("./BZIntegration.jl")
end
