abstract type DummyGrid <: Dispersions.KGridType end
struct DummyKGrid <: FullKGrid{DummyGrid} end
dG = DummyKGrid()

@test_throws MethodError Dispersions.reduceKGrid(dG)
@test_throws ArgumentError gridshape(dG)
