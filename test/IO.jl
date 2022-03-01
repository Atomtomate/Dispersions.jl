@testset "type pretty prints" begin
    io = IOBuffer();
    grid = gen_kGrid("2Dsc-1.0",2)
    print(io, grid)
    @test String(take!(io)) == "SC(t=1.0) grid in 2 dimensions with 4 k-points."
end
