import Base.collect

# ================================================================================ #
#                                Simple Cubic 2D                                    #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type FileDisp <: KGridType end

"""
    FullKGrid_File  <: FullKGrid{FileDisp}

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_File  <: FullKGrid{FileDisp}
    Nk::Int
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
    shape::Tuple{Int,Int,Int}
    function FullKGrid_File(path::String)
        !isfile(path) && error("Dispersion file not found!")
        shape, kGrid, disp = h5open(path, "r") do f
            read(f, "shape"), read(f, "data")[:,1:3], read(f, "data")[:,4]
        end
        new(length(disp), kGrid, disp, shape)
    end
end



"""
    ReducedKGrid_File  <: ReducedKGrid{FileDisp}

For now, no option to read reduction operations from files is implemented.
Therefore the reduced and full k-grids are identical.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_File  <: ReducedKGrid{FileDisp}
    Nk::Int
    kMult::Array{Float64,1}
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
    shape::Tuple{Int,Int,Int}
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::FullKGrid_File) = kG.shape
gridshape(kG::ReducedKGrid_File) = kG.shape

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{FileDisp})

Returns the grid on the fully irredrucible BZ.
For now, reading of symetry operations is not suported and this defaults
to the identity operation.
"""
function reduceKGrid(kG::FullKGrid{FileDisp})
    kmult = ones(kG.Nk)
    return ReducedKGrid_File(kG.Nk, kmult, kG.kGrid, kG.ϵkGrid, kG.shape)
end

expandKArr(kG::ReducedKGrid{FileDisp}, arr::Array{T, 1}) where T = reshape(arr, kG.shape)

reduceKArr(kG::ReducedKGrid{FileDisp}, arr::AbstractArray) = arr[:]

"""
	kGrid_multiplicity_File(kIndices)

Given a set of reduced indices, produce list of multiplicities for each point.
Reading of symmetry operations not yet supported. This will return `ones(length(kIndices))`.
"""
function kGrid_multiplicity_File(kIndices)
    return ones(length(kIndices))
end
