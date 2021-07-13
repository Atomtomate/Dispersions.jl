function gen_kGrid(kg::String, Nk::Int=0; full=false)
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    grid = if lowercase(data[1]) == "3dsc"
        FullKGrid_cP_3D(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "2dsc"
        FullKGrid_cP_2D(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "p6m"
        FullKGrid_p6m(Nk, parse(Float64, data[2]))
    elseif lowercase(data[1]) == "file"
        FullKGrid_File(data[2])
    else
        throw(ArgumentError("Unkown grid type"))
    end
    res = full ? grid : reduceKGrid(grid)
    return res
end

# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #

# ------------------------------ Helper Functions -----------------------------

"""
    gridshape(kG::T) where T <: KGrid

shape of kGrid (e.g. `(kG.Ns, kG.Ns)` for 2D sc) 
"""
gridshape(kG::T) where T <: KGrid = throw(ArgumentError("Unkown kGrid"))

"""
    Nk(kG::T) where T <: KGrid

Total number of k points (length of `kGrid.kGrid` for full grids). 
"""
Nk(kG::T) where T <: KGrid = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

Number of grid points as integer value. `Nk` needs to be accessible for all implemented 
full and reduced k grids.
"""
gridPoints(kG::T) where T <: KGrid = kG.kGrid

"""
    dispersion(kG::T)::Int where T <: KGrid

Returns dispersion relation of grid.
"""
dispersion(kG::T) where T <: KGrid = kG.ϵkGrid

# ------------------------------ Helper Functions -----------------------------
"""
    reduceKGrid(kGrid::FullKGrid{T}) where T <: KGridType

Returns the grid on the fully irredrucible BZ.
"""
reduceKGrid(kGrid::FullKGrid{T}) where T <: KGridType = throw(MethodError("KGrid Instance not found")) 

expandKArr(kGrid::ReducedKGrid{T}, arr::Array) where T <: KGridType = throw(MethodError("KGrid Instance not found")) 


# ------------------------------ Private Functions -----------------------------
"""
    gen_ϵkGrid(::T1, kGrid::T2, t::T3) where {T1 <: KGridType, T2 <: Array, T3 <: Number}

Generates dispersion relation for array of k-points. Usually only accessed by constructor of
KGrid
"""
gen_ϵkGrid(::Type{T1}, kGrid::T2, t::T3) where {T1 <: KGridType, T2 <: Array, T3 <: Number} = throw(MethodError("KGrid Instance not found"))


function conv(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    Nk(kG) == 1 && return arr1 .* arr2
    reshape(fft(expandKArr(kG, arr1)) .* fft(expandKArr(kG, arr2)), gridshape(kG)) |> ifft |> reverse |> x-> reduceKArr(kG, x) ./ Nk(kG)
end


function conv_fft1(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    Nk(kG) == 1 && return (arr1 .* arr2
    reshape(fft(expandKArr(kG, arr1))[:] .* arr2, gridshape(kG)) |> ifft |> reverse |> x-> reduceKArr(kG, x) ./ Nk(kG)
end

function conv_fft(kG::ReducedKGrid, arr1::AbstractArray{Complex{Float64},1}, arr2::AbstractArray{Complex{Float64},1})
    Nk(kG) == 1 && return arr1 .* arr2
    reshape(arr1 .* arr2, gridshape(kG)...) |> ifft |> reverse |> x-> reduceKArr(kG, x) ./ Nk(kG)
end

