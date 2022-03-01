using Test

using Dispersions


grid_list = ["2Dsc-1.3", "3Dsc-1.3", "p6m-1.1"]

include("./helper_functions.jl")
@testset "IO" begin
    include("./Types.jl")
    include("./IO.jl")
end

@testset "KGrid" begin
    include("./KGrid.jl")
end


# @testset "BZ Integration" begin
#     include("./BZIntegration.jl")
# end

# @testset "common" begin
#     include("./common.jl")
# end
