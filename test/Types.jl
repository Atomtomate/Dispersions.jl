import Dispersions.KGridType
abstract type DummyGrid <: KGridType end

#@test_throws MethodError Dispersions.reduceKGrid(dG)
@test_throws ArgumentError Dispersions.gen_sampling(DummyGrid, 0, 0)
@test_throws ArgumentError Dispersions.basis_transform(DummyGrid, [])
@test_throws ArgumentError Dispersions.reduce_KGrid(DummyGrid, 0, 0, [])
#@test_throws ArgumentError Dispersions.gen_ϵkGrid(DummyGrid, [(0,)], 0.0)
