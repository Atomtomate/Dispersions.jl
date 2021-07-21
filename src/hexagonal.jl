import Base.collect

# ================================================================================ #
#                                Simple Cubic 2D                                    #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type p6m <: KGridType end

"""
    FullKGrid_p6m  <: FullKGrid{p6m}

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_p6m  <: FullKGrid{p6m}
    Nk::Int
    Ns::Int
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
    t::Float64
    function FullKGrid_p6m(Nk::Int, t::Float64)
        ki_monkhorst = [(2*j - Nk - 1)/(2*Nk) for j in 1:Nk]
        ki_mesh = Base.product(ki_monkhorst,ki_monkhorst);
        #f_old(x) = (2π * (x[1] - x[2]/sqrt(3)), 4π*x[2]/sqrt(3))
        f(x) =  (2π * x[1], -2π * x[1]/sqrt(3) + 4π*x[2]/sqrt(3))
        kGrid = collect(map(f, ki_mesh))[:]
        disp = gen_ϵkGrid(p6m, kGrid, t)
        new(Nk^2, Nk, kGrid, disp, t)
    end
end

"""
    ReducedKGrid_p6m  <: ReducedKGrid{p6m}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_p6m  <: ReducedKGrid{p6m}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::GridInd2D
    kMult::Array{Float64,1}
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::FullKGrid_p6m) = (kG.Ns, kG.Ns)
gridshape(kG::ReducedKGrid_p6m) = (kG.Ns, kG.Ns)

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{p6m})

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{p6m})
    kGrid = reshape(kG.kGrid, gridshape(kG))
    ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:2]...))

    nh = floor(Int,gridshape(kG)[2]/2) + 1
    index = [[x, y] for y in nh:gridshape(kG)[2] for x in 1:gridshape(kG)[1]]


    ind_red = GridInd2D(undef, length(index))
    grid_red = GridPoints2D(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = kGrid[ti...]
        ϵk_red[i] = ϵkGrid[ti...]
    end
	kmult = kGrid_multiplicity(p6m, ind_red)
    return ReducedKGrid_p6m(kG.Nk, kG.Ns, kG.t, ind_red, kmult, grid_red,  ϵk_red)
end

"""
    expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2,1}) where {T1 <: Union{cP_2D,cP_3D}, T2 <: Any

Expands array of values on reduced k grid back to full BZ.
"""
function expandKArr(kG::ReducedKGrid{p6m}, arr::Array{T, 1}) where T
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    return arr
end

function reduceKArr(kG::ReducedKGrid{p6m}, arr::AbstractArray)
    return arr[:]
end

"""
	kGrid_multiplicity_p6m(kIndices)

Given a set of reduced indices, produce list of multiplicities for each point
"""
function kGrid_multiplicity(::Type{p6m}, kIndices)
    res = ones(Float64, length(kIndices))
    nh = minimum(last.(kIndices))
    for i in axes(res,1)
        (kIndices[i][2] != nh) && (res[i] = 2.0)
    end
    return res
end

gen_ϵkGrid(::Type{p6m}, kGrid::GridPoints2D, t::T1) where T1 <: Number = collect(map(kᵢ -> -2*t*(cos.(0.5*(kᵢ[1] + sqrt(3)*kᵢ[2])) + cos(0.5*(kᵢ[1] - sqrt(3)*kᵢ[2])) + cos(kᵢ[1])), kGrid))
