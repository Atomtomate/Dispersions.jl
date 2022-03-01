abstract type DummyGrid <: Dispersions.KGridType end
struct DummyKGrid end
dG = DummyKGrid()

@test_throws MethodError Dispersions.KGrid(dG)
@test_throws MethodError gridshape(dG)
