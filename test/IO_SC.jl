@testset "type pretty prints" begin
    io = IOBuffer();
    grid = Dispersions.gen_cP_kGrid(4,2, 1.0)
    print(io, ind)
    @test String(take!(io)) == "Index[(1, 1):(4, 4)]: length:16"
    gp = gridPoints(grid)
    print(io, gp)
    @test String(take!(io)) == "GridPoints[(-1.5707963267948966, -1.5707963267948966):(3.141592653589793, 3.141592653589793)]: length:16"
    print(io, grid)
    @test String(take!(io)) == "FullKGrid_cP[4] for 2D"
end
