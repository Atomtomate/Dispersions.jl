import Base.collect
#TODO update tests and docu for improved expandKArr


"""
    gen_cP_kGrid(Nk::Int64, D::Int64, t::Float64)

Generates a simple cubic lattice in `D` Dimensions
See also [`FullKGrid_cP_3D`](@ref) and [`FullKGrid_cP_3D`](@ref)
# Examples
```
julia> gen_cP_kGrid(2, 2)
"""
function gen_cP_kGrid(Nk::Int64, D::Int64, t::Float64)
    if D == 2
        return FullKGrid_cP_2D(Nk, t)
    elseif D == 3
        return FullKGrid_cP_3D(Nk, t)
    else
        throw("Simple Cubic only implemented for 2D and 3D")
    end
end


# ================================================================================ #
#                                Simple Cubic 2D                                    #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type cP_2D <: KGridType end

"""
    FullKGrid_cP_2D  <: FullKGrid{cP_2D}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_cP_2D  <: FullKGrid{cP_2D}
    Nk::Int
    Ns::Int
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
    t::Float64
    fftw_plan::FFTW.cFFTWPlan
    function FullKGrid_cP_2D(Nk::Int, t::Float64; fftw_plan=nothing)
        kx = [(2*π/Nk) * j - π for j in 1:Nk]
        kGrid  = collect(Base.product([kx for Di in 1:2]...))[:]
        fftw_plan = fftw_plan === nothing ? plan_fft!(randn(Complex{Float64}, Nk, Nk), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new(Nk^2, Nk, kGrid, gen_ϵkGrid(cP_2D,kGrid,t),t,fftw_plan)
    end
end

"""
    ReducedKGrid_cP_2D  <: ReducedKGrid{cP_2D}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_cP_2D  <: ReducedKGrid{cP_2D}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::GridInd2D
    kMult::Array{Float64,1}
    kGrid::GridPoints2D
    ϵkGrid::GridDisp
    expand_perms::Vector{Vector{CartesianIndex{2}}}
    expand_cache::Array{Complex{Float64}}
    fftw_plan::FFTW.cFFTWPlan
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::ReducedKGrid_cP_2D) = (kG.Ns, kG.Ns)
gridshape(kG::FullKGrid_cP_2D) = (kG.Ns, kG.Ns)

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{T}) where T <: Union{cP_2D, cP_3D}

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{cP_2D})
    D = length(gridshape(kG))
    kGrid = reshape(kG.kGrid, gridshape(kG))
    ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:D]...))

    ll = floor(Int,size(kGrid,1)/2 + 1)
    la = ceil(Int,kG.Ns/2) .- 1
    index = [[x + la, y + la]  for x=1:ll for y=1:x]

    ind_red = GridInd2D(undef, length(index))
    grid_red = GridPoints2D(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))

    # Compute reduced arrays
    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = kGrid[ti...]
        ϵk_red[i] = ϵkGrid[ti...]
    end

    kMult, expand_perms = build_expand_mapping_SC(D, kG.Ns, ind_red)
    expand_cache = Array{Complex{Float64}}(undef, gridshape(kG))

    return ReducedKGrid_cP_2D(kG.Nk, kG.Ns, kG.t, ind_red, kMult, grid_red, ϵk_red,
                              expand_perms, expand_cache, kG.fftw_plan)
end

"""
    expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2,1}) where {T1 <: Union{cP_2D,cP_3D}, T2 <: Any

Expands array of values on reduced k grid back to full BZ.
"""

function expandKArr(kG::ReducedKGrid{cP_2D}, arr::Array{T, 1}) where T
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    newArr = Array{eltype(arr)}(undef, gridshape(kG)...)
    for (ri,perms) in enumerate(kG.expand_perms)
        for p in perms
            newArr[p] = arr[ri]
        end
    end
    #kG.Ns > 2 && expand_mirror!(newArr)
    return newArr
end

function expandKArr!(kG::ReducedKGrid{cP_2D}, arr::Array{Complex{Float64}, 1})
    for (ri,perms) in enumerate(kG.expand_perms)
        for p in perms
            kG.expand_cache[p] = arr[ri]
        end
    end
    #kG.Ns > 2 && expand_mirror!(kG.expand_cache)
end

# ================================================================================ #
#                                Simple Cubic 3D                                   #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #
abstract type cP_3D <: KGridType end

"""
    FullKGrid_cP_3D  <: FullKGrid{cP_3D}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
"""
struct FullKGrid_cP_3D  <: FullKGrid{cP_3D}
    Nk::Int
    Ns::Int
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
    t::Float64
    fftw_plan::FFTW.cFFTWPlan
    function FullKGrid_cP_3D(Nk::Int, t::Float64; fftw_plan=nothing)
        fftw_plan = fftw_plan === nothing ? plan_fft!(randn(Complex{Float64}, Nk, Nk, Nk), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        kx = [(2*π/Nk) * j - π for j in 1:Nk];
        kGrid  = collect(Base.product([kx for Di in 1:3]...))[:]
        new(Nk^3, Nk, kGrid, gen_ϵkGrid(cP_3D,kGrid,t),t, fftw_plan)
    end
end


"""
    ReducedKGrid_cP_3D  <: ReducedKGrid{cP_3D}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Nk`**    : `Int` k points in full grid.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_cP_3D  <: ReducedKGrid{cP_3D}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::GridInd3D
    kMult::Array{Float64,1}
    kGrid::GridPoints3D
    ϵkGrid::GridDisp
    expand_perms::Vector{Vector{CartesianIndex{3}}}
    expand_cache::Array{Complex{Float64}}
    fftw_plan::FFTW.cFFTWPlan
end

gridshape(kG::ReducedKGrid_cP_3D) = (kG.Ns, kG.Ns, kG.Ns)
gridshape(kG::FullKGrid_cP_3D) = (kG.Ns, kG.Ns, kG.Ns)

function reduceKGrid(kG::FullKGrid{cP_3D})
    kGrid = reshape(kG.kGrid, gridshape(kG))
    ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    D = length(gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:D]...))

    ll = floor(Int,size(kGrid,1)/2 + 1)
    la = ceil(Int,kG.Ns/2) .- 1
    index = [[x+la,y+la,z+la] for x=1:ll for y=1:x for z = 1:y]
    ind_red = GridInd3D(undef, length(index))
    grid_red = GridPoints3D(undef, length(index))
    ϵk_red = GridDisp(undef, length(index))


    for (i,ti) in enumerate(index)
        ind_red[i] = ind[ti...]
        grid_red[i] = kGrid[ti...]
        ϵk_red[i] = ϵkGrid[ti...]
    end

    kMult, expand_perms = build_expand_mapping_SC(D, kG.Ns, ind_red)
    expand_cache = Array{Complex{Float64}}(undef, gridshape(kG))

    return ReducedKGrid_cP_3D(kG.Nk, kG.Ns, kG.t, ind_red, kMult, grid_red, ϵk_red, 
                              expand_perms, expand_cache, kG.fftw_plan)
end


function expandKArr(kG::ReducedKGrid{T1}, arr::Array{T2, 1}) where {T1 <: Union{cP_2D, cP_3D}, T2 <: Any}
    length(arr) != length(kG.kInd) && throw(ArgumentError("length of k grid ($(length(kG.kInd))) and argument ($(length(arr))) not matching"))
    newArr = Array{T2}(undef, gridshape(kG))
    for (ri,perms) in enumerate(kG.expand_perms)
        for p in perms
            newArr[p] = arr[ri]
        end
    end
    return newArr
end

function expandKArr!(kG::ReducedKGrid{T1}, arr::Array{Complex{Float64}, 1}) where {T1 <: Union{cP_2D, cP_3D}}
    for (ri,perms) in enumerate(kG.expand_perms)
        for p in perms
            kG.expand_cache[p] = arr[ri]
        end
    end
end


function reduceKArr(kG::ReducedKGrid{T1}, arr::AbstractArray) where {T1 <: Union{cP_2D, cP_3D}}
    res = Array{eltype(arr), 1}(undef, length(kG.kInd))
    for (i,ki) in enumerate(kG.kInd)
        res[i] = arr[ki...]
    end
    return res
end

function reduceKArr!(kG::ReducedKGrid{T1}, res::AbstractArray, arr::AbstractArray) where {T1 <: Union{cP_2D, cP_3D}}
    for (i,ki) in enumerate(kG.kInd)
        res[i] = arr[ki...]
    end
end


#TODO: this is a placeholder until convolution is ported from lDGA code
function reduceKArr_reverse(kG::ReducedKGrid{T1}, arr::AbstractArray) where {T1 <: Union{cP_2D, cP_3D}}
    res =  Array{eltype(arr)}(undef, length(kG.kInd))
    N = size(arr)
    ll = floor(Int,size(arr,1)/2 + 1)
    index = T1 === cP_2D ? [[x,y] for x=1:ll for y=1:x] : [[x,y,z] for x=1:ll for y=1:x for z = 1:y]
    for (i,ti) in enumerate(index)
        res[i] = arr[(N .- ti)...]
    end
    return res
end


gen_ϵkGrid(::Type{cP_2D}, kGrid::GridPoints2D, t::T1) where T1 <: Number = collect(map(kᵢ -> -2*t*sum(cos.(kᵢ)), kGrid))
gen_ϵkGrid(::Type{cP_3D}, kGrid::GridPoints3D, t::T1) where T1 <: Number = collect(map(kᵢ -> -t*sum(cos.(kᵢ)), kGrid))
ifft_post!(::Type{cP_2D}, x::Array{T,2}) where T <: Number = reverse!(x) 
ifft_post!(::Type{cP_3D}, x::Array{T,3}) where T <: Number = reverse!(x) 

function build_expand_mapping_SC(D::Int, Ns::Int, ind_red::Array)
    expand_perms = Vector{Vector{CartesianIndex{D}}}(undef, length(ind_red))
    kMult = Array{Int,1}(undef, length(ind_red))

    mirror_list = Array{NTuple{D,Int},1}()
    for i in 1:D
        push!(mirror_list, map( x-> tuple((x)...) ,unique(permutations([(j <= i) ? Ns : 0 for j in 1:D ])))...)
    end
    #  - Expand mapping
    for (ri, redInd) in enumerate(ind_red)
        perms = unique(permutations(redInd))
        expand_perms[ri] = Vector{CartesianIndex{D}}()
        for (ip,p) in enumerate(perms)
            push!(expand_perms[ri], CartesianIndex(p...))
            for mi in mirror_list
                if all(abs.(mi .- p) .> 0)
                    push!(expand_perms[ri], CartesianIndex(abs.(mi .- p)...))
                end
            end
        end
        expand_perms[ri] = unique(expand_perms[ri])
        kMult[ri] = length(expand_perms[ri])
    end
    return kMult, expand_perms
end
