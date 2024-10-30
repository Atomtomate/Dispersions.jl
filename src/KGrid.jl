#TODO: t/tp/tpp should be a vector
# TODO: rename gen_shifted_ϵkGrid to something appropriate

# ============================================ KGrid Type ============================================
"""
    KGrid{T <: KGridType, D}

Fields
-------------
- **`Nk`**      : `Int`, Number of total k-points
- **`Ns`**      : `Int`, Number of sampling points per dimension
- **`t`**       : `Float64`, hopping parameter
- **`tp`**      : `Float64`, hopping parameter, next nearest neighbor
- **`tpp`**     : `Float64`, hopping parameter, next next nearest neighbor
- **`kGrid`**   : `Vector{NTuple{D,Float64}}`, vector of k-points. Each element is a D-tuple
- **`ϵkGrid`**  : `Vector{Float64}`, Dispersion relation
- **`kInd`**    : `Vector{NTuple{D,Int}}`, vector of indices mapping from the full to reduced lattice.
- **`kInd_conv`**    : `Vector{NTuple{D,Int}}`, vector of indices mapping from the full to reduced lattice after convolution (this incorporates possible reorderings).
- **`kInd_crossc`**   : `Vector{NTuple{D,Int}}`, vector of indices mapping from the full to reduced lattice after crosscorrelation (this incorporates possible reorderings).
- **`kMult`**   : `Vector{Int}`, multiplicity per k-point in reduced lattice
- **`expand_perms`** : `Vector{NTuple{D, Int}}`, mapping of each k-point in reduced lattice to full lattice points
- **`expand_cache`** : `Array{ComplexF64}`, internal cache for expansion of reduced to full lattice before executing convolutions
- **`conv_cache`**   : `Array{ComplexF64,D}`, innternal cache for convolutions
- **`fftw_plan`**    : `FFTW.cFFTWPlan`, fft plan to be executed in convolutions. WARNING: This field can not be serialized right now and needs to be reconstructed after reading a `KGrid` from disk.
"""
struct KGrid{T <: KGridType, D}
    Nk::Int
    Ns::Int
    t::Float64
    tp::Float64
    tpp::Float64
    k0::Tuple
    kGrid::GridPoints
    ϵkGrid::GridDisp
    kInd::GridInd
    kInd_conv::GridInd
    kInd_crossc::GridInd
    kMult::Array{Float64,1}
    expand_perms::Vector{Vector{CartesianIndex{D}}}
    cache1::Array{ComplexF64,D}
    cache2::Array{ComplexF64,D}
    fftw_plan::FFTW.cFFTWPlan
    function KGrid(GT::Type{T}, D::Int, Ns::Int, t::Float64, tp::Float64, tpp::Float64; full::Bool=false, fftw_plan=nothing) where T<:KGridType
        sampling = gen_sampling(GT, D, Ns)
        kGrid_f = map(v -> basis_transform(GT, v), sampling)
        kInd, kInd_conv, kInd_crossc, kMult, expand_perms, kGrid = full ? reduce_KGrid_ident(GT, D, Ns, kGrid_f) : reduce_KGrid(GT, D, Ns, kGrid_f)
        ϵkGrid =  gen_ϵkGrid(GT, kGrid, t, tp, tpp)
        gs = repeat([Ns], D)
        k0 = findfirst(k -> all(isapprox.(k, 0, atol=1e-7)), kGrid_f)
        if isnothing(k0) 
            #println("WARNING!!! k-grid sampling must contain zero-vector in order for convolutions to work!")
            k0 = first(CartesianIndices([1 2; 3 4]))
        end
        fftw_plan = fftw_plan === nothing ? plan_fft!(FFTW.FakeArray{ComplexF64}(gs...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{GT,D}(Ns^D, Ns, t, tp, tpp, k0, kGrid, ϵkGrid, kInd, kInd_conv, kInd_crossc, kMult, expand_perms,
                  Array{ComplexF64,D}(undef, gs...), Array{ComplexF64,D}(undef, gs...), fftw_plan)
    end
end

"""
    gen_kGrid(kg::String, Ns::Int; full=false)

Generates a KGrid of type and hopping strength, given in `kG` with `Ns` sampling points in the first Brillouin zone. Options are:
- '3dcP-...'         : simple cubic 3D
- '2dcP-...'         : simple cubic 2D
- '2dcP-...-...-...' : simple cubic 2D with next-next nearest neighbor hopping
- 'cF-...'           : FCC
- 'cI-...'           : BCC

If `full=true` is set, no symmetry reduction is used.

# Examples
```
julia> gen_kGrid("3dcP-1.5", 10)
cP(t=1.5) grid in 3 dimensions with 1000 k-points.
```
"""
function gen_kGrid(kg::String, Ns::Int; full=false)
    findfirst("-", kg) === nothing && throw(ArgumentError("Please provide lattice type and hopping, e.g. SC3D-1.1"))
    tp = 0.0
    tpp = 0.0
    str_v = split(kg,"-")
    data = deepcopy(str_v)
    ii = findall(isempty, str_v)
    for i in ii
        data[i+1] = string("-",data[i+1])
    end
    deleteat!(data, ii)
    gt_s = lowercase(data[1])
    gt_s = replace(gt_s, "cp" => "sc", "cf" => "fcc", "ci" => "bcc")
    t = parse(Float64, data[2])
    if length(data) == 3
        tp = parse(Float64, data[3])
        gt_s = endswith(gt_s, "sc") ? string(gt_s, "nn") : gt_s
    elseif  length(data) == 4
        tp = parse(Float64, data[3])
        tpp = parse(Float64, data[4])
        gt_s = endswith(gt_s, "sc") ? string(gt_s, "nn") : gt_s
    end
    if gt_s == "3dsc"
        KGrid(cP, 3, Ns, t, tp, tpp; full=full)
    elseif gt_s == "2dscnn"
        KGrid(cPnn, 2, Ns, t, tp, tpp; full=full)
    elseif gt_s == "2dsc"
        KGrid(cP, 2, Ns, t, tp, tpp; full=full)
    elseif gt_s == "4dsc" || gt_s == "4dscnn"
        KGrid(cP, 4, Ns, t, tp, tpp; full=full)
    elseif gt_s == "fcc"
        KGrid(cF, 3, Ns, t, tp, tpp; full=full)
    elseif gt_s == "bcc"
        KGrid(cI, 3, Ns, t, tp, tpp; full=full)
    elseif startswith(gt_s, "hofstadter")
        println(gt_s, ", ", t, ", ", tp, ", ", tpp)
        GridName,P,Q = split(gt_s, ":")
        P = parse(Int, P)
        Q = parse(Int, Q)
        KGrid(Hofstadter{P,Q}, 2, Ns, t, tp, tpp)
    else
        throw(ArgumentError("Unkown grid type: $kg / $gt_s"))
    end
end


"""
    transform_to_first_BZ(kG, k)

Returns a new vector `kp` that is transformed back into the sample region of `kG`.
"""
transform_to_first_BZ(kG::KGrid{GridType,D}, k) where {GridType,D} = 
    throw(ArgumentError("Cannot transform back to first BZ! Grid type $(typeof(kG)) unkown!"))


"""
    ϵ_k_plus_q(kG::KGrid, q::NTuple)
    
    Evaluates the dispersion relation on the given reciprocal space but expanded and shifted by a constant vector `q`. The corresponding points in reciprocal space are given by `expandKArr(kG, gridPoints(kG))`.

Returns:
-------------
``\\epsilon(k+\\mathrm{shift}): `Vector{NTuple{D,Float64}}`, where D is the diemenion of the grid. Dispersion relation evaluated on the given grid but shifted by the the vector q.

ATTENTION: So far this function is tested for the simple cubic lattice only!

Arguments:
-------------
- `kG`       : reciprocal lattice
- **`q`**    : vector in reciprocal space
"""
function ϵ_k_plus_q(kG::KGrid, q::NTuple)
    if grid_dimension(kG) != length(q)
        throw(ArgumentError("Grid dimension differs from shift dimension!"))
    else
        k_plus_q = map(ki -> ki .+ q, expandKArr(kG, gridPoints(kG))[:])
        return gen_ϵkGrid(grid_type(kG), k_plus_q, kG.t)
    end
end

"""
    grid_type(kG::KGrid)

Returns the grid type of `kG`.

Used to call functions like [`gen_sampling`](@ref gen_sampling), that expect a grid type as input.
"""
function grid_type(kG::KGrid)
    return typeof(kG).parameters[1]
end

"""
    grid_dimension(kG::KGrid)

Returns the grid dimension of `kG`.

Used to call functions like [`gen_sampling`](@ref gen_sampling), that expect a grid dimension as input.
"""
function grid_dimension(kG::KGrid)
    return typeof(kG).parameters[2]
end

"""
    map_to_indices(path::AbstractVector, grid::AbstractArray)

Maps vectors in the given input array `path` to indices in the given k-grid `grid`. 
This is used internally, to construct indices for symmetry paths.
"""
function map_to_indices(path::AbstractVector, grid::AbstractArray)
    result_points = Array{CartesianIndex,1}(undef, length(path))
    residual_vals = Array{Float64,1}(undef, length(path))
    for (i,point) in enumerate(path)
        residual_vals[i], result_points[i] = findnearest(point, grid)
    end
    return result_points, residual_vals
end

"""
    reduce_KGrid_ident(Ns::Int, kGrid::AbstractArray)

Identity operation for KGrid reduction. 
    Used if a full KGrid is requested by providing [`gen_kGrid`](@ref gen_kGrid) with `full=true`.
"""
function reduce_KGrid_ident(gt::GT, D::Int, Ns::Int, kGrid::AbstractArray) where GT
    ind = collect(Base.product([1:Ns for Di in 1:D]...))
    index = [CartesianIndex(el) for el in ind[:]]
    kMult = ones(Int, length(ind))
    expand_perms = map(x -> [CartesianIndex{D}(x)],ind[:])
    k0 = floor.(Int, Tuple(repeat([Ns],D)) ./ 2) .- 1
    m1 = -1 .* k0 .+ 0
    k0, m1 = conv_Indices(gt, D, Ns)
    ind_red_conv  = CartesianIndex.(circshift(ind, m1)[index]); # indices after conv
    ind_red_crossc = CartesianIndex.(circshift(reverse(ind), k0)[index]); # indices after crossc
    # Change from CartesianIndices to LinearIndices for performance reasons
    I = LinearIndices(ind)
    index = I[index]
    ind_red_conv = I[ind_red_conv]
    ind_red_crossc = I[ind_red_crossc]
    return index, ind_red_conv, ind_red_crossc, kMult, expand_perms, kGrid[:]
end

# ================================= Subsampling (Reimplement this!!!) ================================
"""
    check_subsample(sub_sample, full_sample)

Internal function, checking if all sample points of a sub-sample are contained in full sample.
See [`build_kGrid_subsample`](@ref build_kGrid_subsample).
"""
function check_subsample(sub_sample, full_sample)
    test_k_subsample = true
    for kv in sub_sample
        if !(kv in full_sample)
            test_k_subsample = false
        end
    end
    return test_k_subsample
end

subsample_indices(sub_sample, full_sample) = map(kv -> findfirst(el -> all(el .≈ kv),full_sample), sub_sample)

"""
    build_kGrid_subsample(kG_full::Kgrid, Nk_sub::Int)

Builds new [`KGrid`](@ref KGrid) of the same type as `kG_full` with at most `Nk_sub` sample points.

Returns new `KGrid` and `Vector{Cartesianindex}` of indices for indices corresponding to full grid.

TODO: this is brute-force checking all sizes for now!
"""
function build_kGrid_subsample(kG_full::KGrid, Nk_sub::Int)
    println("Brute force checking for possible sub-samples. This needs to be rewritten")
    full_sample_approx = map(x -> round.(x, digits=5), kG_full.kGrid)
    res = nothing
    sub_indices = nothing
    
    for Nk_test = Nk_sub:-1:1
        ti = if iseven(Nk_test)
            kG_sub = KGrid(grid_type(kG_full), grid_dimension(kG_full), Nk_test, kG_full.t, kG_full.tp, kG_full.tpp)
            sub_sample_approx = map(x -> round.(x, digits=5), kG_sub.kGrid)
            ti = check_subsample(sub_sample_approx, full_sample_approx)
            if ti 
                res = kG_sub
                sub_indices = subsample_indices(sub_sample_approx, full_sample_approx)
            end
            ti
        else
            false
        end 
        ti && break
    end
    return res, sub_indices
end
