using Test

using Dispersions

include("./helper_functions.jl")
@testset "SC" begin
    include("./Types.jl")
    include("./SC.jl")    
    #include("./IO_SC.jl")    
end

@testset "BZ Integration" begin
    include("./BZIntegration.jl")
end
