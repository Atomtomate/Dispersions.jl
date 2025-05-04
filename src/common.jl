# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #
export reduceKArr,
    reduceKArr!,
    expandKArr,
    expandKArr!,
    reverseKArr,
    conv,
    conv!,
    conv_fft,
    conv_fft!,
    conv_fft1,
    conv_fft1!

# ------------------------------ Access Functions -----------------------------
"""
    gridshape(kG::T) where T <: KGrid

shape of kGrid (e.g. `(kG.Ns, kG.Ns)` for 2D cP) 
"""
gridshape(kG::KGrid{T,D}) where {T,D} = ntuple(_ -> kG.Ns, D)

"""
    Nk(kG::T) where T <: KGrid

Total number of k points (length of `kGrid.kGrid` for full grids). 
"""
Nk(kG::KGrid)::Int = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

k vektors for the given grid. Elements of the irreducible part only.
"""
gridPoints(kG::KGrid) = kG.kGrid

"""
    dispersion(kG::T)::Int where T <: KGrid

Returns dispersion relation of grid.
"""
dispersion(kG::KGrid) = kG.ϵkGrid

# ------------------------- Sampling Related Functions ------------------------
"""
    reduceKArr(kGrid::KGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a full BZ and returns an array with data over fully irreducible
BZ. This is mainly used after convolutions, since they require data over the full BZ.
"""
function reduceKArr(kG::KGrid{gT,D}, arr::AbstractArray{T,D}) where {gT <: KGridType,T,D}
    res = similar(arr, length(kG.kInd))
    reduceKArr!(kG, res, arr)
    return res
end


"""
    reduceKArr!(kGrid::KGrid, res, arr)

Inplace version of [`reduceKArr`](@ref). The results are written to `res`.
"""
function reduceKArr!(
    kG::KGrid{gT,D},
    res::AbstractArray{T,1},
    arr::AbstractArray{T,D},
) where {gT <: KGridType,T,D}
    @inbounds res[:] = arr[kG.kInd]
end


"""
    expandKArr(kGrid::KGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a reduced BZ and returns an array with data over
full BZ. This is mainly used before convolutions, since they require data over the full BZ.
"""
function expandKArr(
    kG::KGrid{gT,D},
    arr::AbstractArray{T,1}; 
)::AbstractArray{T,D} where {gT <: KGridType,T,D}
    length(arr) != length(kG.kInd) && throw(
        ArgumentError(
            "length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching",
        ),
    )
    res = similar(arr, gridshape(kG)...)
    expandKArr!(kG, res, arr)
    return res
end

"""
    expandKArr!(kG, arr)

Inplace version of [`expandKArr`](@ref). The results are by default written to `kG.cache1`.
"""
function expandKArr!(
    kG::KGrid, 
    arr::Array{ComplexF64, 1}; 
    res::Array{ComplexF64, D} = kG.cache1
)::Nothing where {D}
    expandKArr!(kG, res, arr)
    return nothing
end

"""
    expandKArr!(kG, [res,] arr)

Inplace version of [`expandKArr`](@ref). The results are written to `res`.
"""
function expandKArr!(
    kG::KGrid{gT,D},
    res::AbstractArray{T,D},
    arr::Array{T,1},
)::Nothing where {gT <: KGridType,T,D}
    for (ri, perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds res[p] = arr[ri]
        end
    end
    return nothing
end

"""
    reverseKArr(kGrid::KGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a reduced OR full BZ and returns an array with reversed
k-indices. I.e., ``\\f_q \\to f_{-q}``.
"""
function reverseKArr(
    kG::KGrid{gT,D},
    arr::AbstractArray{T,1},
)::AbstractArray{T,1} where {gT <: KGridType,T,D}
    reverseKArr(kG, expandKArr(kG, arr))
end

function reverseKArr(
    kG::KGrid{gT,D},
    arr::AbstractArray{T,D},
)::AbstractArray{T,D} where {gT <: KGridType,T,D}
    shift_vec = 2 .* kG.k0 .- gridshape(kG) .- 1
    return circshift(reverse(arr), shift_vec)
end

# ------------------------------ Convolution Functions -----------------------------
"""
    conv(kG::KGrid, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Computes the convolution with a plus sign of data over two arrays `arr1` and `arr2`,
i.e. ``res[k] = \\sum_{q \\in \\text{BZ}} arr1[k] * arr2[k+q]``.

`crosscorrelation` sets the sign of the 'convolution' to `+`, i.e. ``\\sum_j f_i g_{i+j}`` instead of ``\\sum_j f_i g_{i-j}``.
"""
function conv(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64,1};
    crosscorrelation::Bool=true,
    cache1::Array{ComplexF64, D} = kG.cache1,
    cache2::Array{ComplexF64, D} = kG.cache2
) where D
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation, cache1=cache1, cache2=cache2)
    return res
end
#TODO: only for testing, remove?
#conv(kG::KGrid, arr1::AbstractArray, arr2::AbstractArray) = conv(kG, convert.(ComplexF64, arr1[:]), convert.(ComplexF64, arr2[:]))

"""
    conv!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv`](@ref). The results are written to `res`.
"""
function conv!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64,1};
    crosscorrelation::Bool=true,
    cache1::Array{ComplexF64, D} = kG.cache1,
    cache2::Array{ComplexF64, D} = kG.cache2
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    expandKArr!(kG, cache2, arr2)
    crosscorrelation && reverse!(cache2)
    kG.fftw_plan * cache2

    conv_fft1!(kG, res, arr1, cache2, crosscorrelation=crosscorrelation, cache=cache1)
end

"""
    conv_fft1(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr2`. Note, that reversing the array is also
expected, i.e. `arr2 = fft(reverse(in_arr2))`, if `crosscorrelation = true` is assumed.
"""
function conv_fft1(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool=true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_fft1!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation, cache=cache)
    return res
end

"""
    conv_fft1!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv_fft1`](@ref).
"""
function conv_fft1!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool=true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    expandKArr!(kG, cache, arr1)
    kG.fftw_plan * cache

    conv_fft!(kG, res, cache, arr2, crosscorrelation=crosscorrelation, cache=cache)
end

"""
    conv_fft(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr1` and `arr2`. Note, that reversing `arr2` is also
expected, see also [`conv_fft1`](@ref), if `crosscorrelation = true` is assumed.
"""
function conv_fft(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool=true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return arr1 .* arr2
    res = Array{ComplexF64,1}(undef, length(kG.kMult))
    conv_fft!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation, cache=cache)
    return res
end


"""
    conv_fft!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv_fft`](@ref).
"""
function conv_fft!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool=true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    for i in eachindex(cache)
        @inbounds cache[i] = arr1[i] * arr2[i]
    end
    kG.fftw_plan \ cache
    conv_post!(kG, res, cache, crosscorrelation=crosscorrelation)
end

# ===== no plan functions. TODO: replace definition with macro + find a way to serialize plans
function conv_noPlan(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64,1};
    crosscorrelation::Bool = true
)
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_noPlan!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation)
    return res
end

function conv_noPlan!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64,1};
    crosscorrelation::Bool = true,
    cache1::Array{ComplexF64, D} = kG.cache1,
    cache2::Array{ComplexF64, D} = kG.cache2
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    expandKArr!(kG, cache2, arr2)
    crosscorrelation && reverse!(cache2)

    fft!(cache2)
    conv_fft1_noPlan!(kG, res, arr1, cache2, crosscorrelation=crosscorrelation, cache=cache1)
end

function conv_fft1_noPlan(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool = true
)
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_fft1!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation)
    return res
end

function conv_fft1_noPlan!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool = true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    expandKArr!(kG, cache, arr1)
    fft!(cache)

    conv_fft_noPlan!(kG, res, cache, arr2, crosscorrelation=crosscorrelation, cache=cache)
end

function conv_fft_noPlan(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool = true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return arr1 .* arr2
    res = Array{ComplexF64,1}(undef, length(kG.kMult))
    conv_fft_noPlan!(kG, res, arr1, arr2, crosscorrelation=crosscorrelation, cache=cache)
    return res
end

function conv_fft_noPlan!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64};
    crosscorrelation::Bool = true,
    cache::Array{ComplexF64, D} = kG.cache1
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    for i in eachindex(cache)
        @inbounds cache[i] = arr1[i] * arr2[i]
    end
    kG.ifftw_plan * cache
    conv_post!(kG, res, cache, crosscorrelation=crosscorrelation)
end

# ------------------------------ Auxiliary Convolution Functions -----------------------------
"""
    conv_post(kG::KGrid, x::Array{T,D})

Post convolution steps e.g. reversing or shifting the result. Some lattice types may overload this, depending
on the sample points. See cP.jl for an example.
`crosscorrelation` sets the `sign` for the convolution. See also  [`conv`](@ref conv)
"""
function conv_post(kG::KGrid, x::Array{T,D}; crosscorrelation=true) where {D,T}
    res = Array{T,1}(undef, length(kG.kMult)) 
    conv_post!(kG, res, x, crosscorrelation=crosscorrelation)
end

"""
    conv_post!(kG::KGrid{cP,D}, res::Array{T,1}, x::Array{T,D}) where {D,T} 

"""
function conv_post!(kG::KGrid, res::AbstractArray{T,1}, x::AbstractArray{T}; crosscorrelation=true) where T 
    norm = Nk(kG)
    res[:] = if crosscorrelation 
        x[kG.kInd_crossc] ./ norm 
    else
        x[kG.kInd_conv] ./ norm
    end
end

"""
    conv_post_add!(kG::KGrid{cP,D}, res::Array{T,1}, x::Array{T,D}) where {D,T} 

Inplace version of [`conv_post`](@ref), but add values instead of replacing them.
"""
function conv_post_add!(kG::KGrid, res::AbstractArray{T,1}, x::AbstractArray{T}; crosscorrelation=true) where T 
    norm = Nk(kG)
    res[:] += crosscorrelation ? x[kG.kInd_crossc] ./ norm : x[kG.kInd_conv] ./ norm
end

# --------------------------- Symmetry Path Sampling --------------------------
"""
    findnearest(p, A::AbstractArray) = findmin(map(vi -> norm(vi .- p),A))

Finds nearest sampling point in grid `A` to point `p`.
"""
findnearest(p::T, A::AbstractArray{T}) where T = findmin(map(vi -> norm(vi .- p),A))

"""
    map_to_indices_full(path::AbstractVector, grid::AbstractArray)

Finds indices for points along path. 
Also returns residual values for points, i.e. norm of distance vector between point on path and point in grid.
"""
function map_to_indices_full(path::AbstractVector, grid::AbstractArray)
    result_points = Array{CartesianIndex,1}(undef, length(path))
    residual_vals = Array{Float64,1}(undef, length(path))
    for (i,point) in enumerate(path)
        residual_vals[i], result_points[i] = findnearest(point, grid)
    end
    return result_points, residual_vals
end

"""
    map_to_indices(path::AbstractVector, grid::AbstractArray, kG::KGrid)

Finds indices for points along path for reduced k-grid. 
Also returns residual values for points, i.e. norm of distance vector between point on path and point in grid.
"""
function map_to_indices(path::AbstractVector, grid::AbstractArray, kG::KGrid)
    result_points, residual_vals = map_to_indices_full(path, grid)
    rr = map(pi -> findfirst(ki_l -> pi in ki_l, kG.expand_perms), result_points)
    return rr, residual_vals
end

"""
    sample_along_path(data::AbstractArray, path::AbstractVector, kG::KGrid)

Sample `data` from reduced grid `kG` along `path`.
Returns data and residual values (`Vector` of distances between points along path and sample points used).
"""
function sample_along_path(data::AbstractArray, path::AbstractVector, kG::KGrid)
    sampling = gen_sampling(grid_type(kG), grid_dimension(kG), kG.Ns)
    grid = map(v -> basis_transform(grid_type(kG), v), sampling);
    result_points, residual_vals = map_to_indices(path, grid, kG)
    return data[result_points], residual_vals
end
