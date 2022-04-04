# ================================================================================ #
#                                   Interface                                      #
# ================================================================================ #
export reduceKArr,
    reduceKArr!,
    expandKArr,
    expandKArr!,
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
Nk(kG) = kG.Nk

"""
    gridPoints(kG::T)::Int where T <: KGrid

Number of grid points as integer value. `Nk` needs to be accessible for all implemented 
full and reduced k grids.
"""
gridPoints(kG::T) where {T<:KGrid} = kG.kGrid

"""
    dispersion(kG::T)::Int where T <: KGrid

Returns dispersion relation of grid.
"""
dispersion(kG::T) where {T<:KGrid} = kG.ϵkGrid

# ------------------------------ Helper Functions -----------------------------
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
    for (i, ki) in enumerate(kG.kInd)
        @inbounds res[i] = arr[ki...]
    end
end


"""
    expandKArr(kGrid::KGrid, arr)

Takes a kGrid `kGrid` and arbitrary data `arr` over a reduced BZ and returns an array with data over
full BZ. This is mainly used before convolutions, since they require data over the full BZ.
"""
function expandKArr(
    kG::KGrid{gT,D},
    arr::AbstractArray{T,1},
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

Inplace version of [`expandKArr`](@ref). The results are written to `kG.expand_cache`.
"""
function expandKArr!(kG::KGrid, arr::Array{Complex{Float64},1})
    for (ri, perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds kG.expand_cache[p] = arr[ri]
        end
    end
end

"""
    expandKArr!(kG, [res,] arr)

Inplace version of [`expandKArr`](@ref). The results are written to `res`.
"""
function expandKArr!(
    kG::KGrid{gT,D},
    res::AbstractArray{T,D},
    arr::Array{T,1},
) where {gT <: KGridType,T,D}
    for (ri, perms) in enumerate(kG.expand_perms)
        @simd for p in perms
            @inbounds res[p] = arr[ri]
        end
    end
end


# ------------------------------ Convolution Functions -----------------------------
"""
    conv(kG::KGrid, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Computes the convolution with a plus sign of data over two arrays `arr1` and `arr2`,
i.e. ``res[k] = \\sum_{q \\in \\text{BZ}} arr1[k] * arr2[k+q]``.
"""
function conv(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64,1},
    arr2::AbstractArray{ComplexF64,1},
)
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv!(kG, res, arr1, arr2)
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
    arr2::AbstractArray{ComplexF64,1},
)
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    gs = gridshape(kG)

    expandKArr!(kG, kG.cache2, arr2)
    reverse!(kG.cache2)
    kG.fftw_plan * kG.cache2

    conv_fft1!(kG, res, arr1, kG.cache2)
end

"""
    conv_fft1(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr2`. Note, that reversing the array is also
expected, i.e. `arr2 = fft(reverse(in_arr2))`.
"""
function conv_fft1(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64},
)
    Nk(kG) == 1 && return arr1 .* arr2
    res = similar(arr2)
    conv_fft1!(kG, res, arr1, arr2)
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
    arr2::AbstractArray{ComplexF64},
)
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)
    expandKArr!(kG, kG.cache1, arr1)
    kG.fftw_plan * kG.cache1

    conv_fft!(kG, res, kG.cache1, arr2)
end

"""
    conv_fft(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Version of [`conv`](@ref) for precomputed `arr1` and `arr2`. Note, that reversing `arr2` is also
expected, see also [`conv_fft1`](@ref).
"""
function conv_fft(
    kG::KGrid,
    arr1::AbstractArray{ComplexF64},
    arr2::AbstractArray{ComplexF64},
)
    Nk(kG) == 1 && return arr1 .* arr2
    res = Array{ComplexF64,1}(undef, length(kG.kMult))
    conv_fft!(kG, res, arr1, arr2)
    return res
end


"""
    conv_fft!(kG::KGrid, res::AbstractVector{ComplexF64}, arr1::AbstractVector{ComplexF64}, arr2::AbstractVector{ComplexF64})

Inplace version of [`conv_fft`](@ref).
"""
function conv_fft!(
    kG::KGrid,
    res::AbstractArray{ComplexF64,1},
    arr1::Array{ComplexF64,D},
    arr2::Array{ComplexF64,D},
) where D
    Nk(kG) == 1 && return (res[:] = arr1 .* arr2)

    for i in eachindex(kG.cache1)
        @inbounds kG.cache1[i] = arr1[i] * arr2[i]
    end
    kG.fftw_plan \ kG.cache1
    conv_post!(kG, res, kG.cache1)
end


# ------------------------------ Auxiliary Convolution Functions -----------------------------
"""
    conv_post(kG::KGrid, x::Array{T,D})

Post convolution steps i.e. reversing the result. Some lattice typs may overload this, depending
on the sample points. See cP.jl for an example.
"""
function conv_post(kG::KGrid, x::Array{T,D}) where {D,T}
    res = Array{T,1}(undef, length(kG.kMult)) 
    conv_post!(kG, res, x)
end

"""
    conv_post!(kG::KGrid{cP,D}, res::Array{T,1}, x::Array{T,D}) where {D,T} 

Inplace version of [`conv_post`](@ref). Warning: `res` should not alias `kG.cache2` as some
implementations may use this cache without explicitly checking for pointer aliases. 
"""
function conv_post!(kG::KGrid, res::AbstractArray{T,1}, x::Array{T,D}) where {D,T} 
    reverse!(x)
    reduceKArr!(kG, res, x)
    norm = Nk(kG)
    res[:] = res ./ norm
end
