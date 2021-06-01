@testset "constructor" begin
    p6m_2d_2 = gen_kGrid("p6m-1.3",2)
    p6m_2d_3 = gen_kGrid("p6m-1.3",4)
    p6m_2d_16 = gen_kGrid("p6m-1.4",16)
    r2 = reduceKGrid(p6m_2d_2)
    r3 = reduceKGrid(p6m_2d_2)
    r16 = reduceKGrid(p6m_2d_16)
    @test Nk(p6m_2d_2) == 2^2
    @test_throws ArgumentError expandKArr(r16, [1,2,3,4])
end

