#TODO: t/tp/tpp should be a vector

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
    kGrid::GridPoints
    ϵkGrid::GridDisp
    kInd::GridInd
    kInd_conv::GridInd
    kMult::Array{Float64,1}
    expand_perms::Vector{Vector{CartesianIndex{D}}}
    cache1::Array{ComplexF64,D}
    cache2::Array{ComplexF64,D}
    fftw_plan::FFTW.cFFTWPlan
    function KGrid(GT::Type{T}, D::Int, Ns::Int, t::Float64, tp::Float64, tpp::Float64; fftw_plan=nothing) where T<:KGridType
        sampling = gen_sampling(GT, D, Ns)
        kGrid_f = map(v -> basis_transform(GT, v), sampling)
        kInd, kInd_conv, kMult, expand_perms, kGrid = reduce_KGrid(GT, D, Ns, kGrid_f)
        ϵkGrid =  gen_ϵkGrid(GT, kGrid, t, tp, tpp)
        gs = repeat([Ns], D)
        fftw_plan = fftw_plan === nothing ? plan_fft!(FFTW.FakeArray{ComplexF64}(gs...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{GT,D}(Ns^D, Ns, t, tp, tpp, kGrid, ϵkGrid, kInd, kInd_conv, kMult, expand_perms,
                  Array{ComplexF64,D}(undef, gs...), Array{ComplexF64,D}(undef, gs...), fftw_plan)
    end
end

"""
    gen_kGrid(kG::String, Ns::Int)

Generates a KGrid of type and hopping strength, given in `kG` with `Ns` sampling points in the first Brillouin zone. Options are:
- '3dcP-...'         : simple cubic 3D
- '2dcP-...'         : simple cubic 2D
- '2dcP-...-...-...' : simple cubic 2D with next-next nearest neighbor hopping
- 'cF-...'           : FCC
- 'p6m-...'          : hexagonal

# Examples
```
julia> gen_kGrid("3dcP-1.5", 10)
cP(t=1.5) grid in 3 dimensions with 1000 k-points.
```
"""
function gen_kGrid(kg::String, Ns::Int)
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
        KGrid(cP, 3, Ns, t, tp, tpp)
    elseif gt_s == "2dscnn"
        KGrid(cPnn, 2, Ns, t, tp, tpp)
    elseif gt_s == "2dsc"
        KGrid(cP, 2, Ns, t, tp, tpp)
    elseif gt_s == "fcc"
        KGrid(cF, 3, Ns, t, tp, tpp)
    elseif gt_s == "bcc"
        KGrid(cI, 3, Ns, t, tp, tpp)
    elseif gt_s == "p6m"
        KGrid(p6m, 2, Ns, t, tp, tpp)
    elseif startswith(gt_s, "hofstadter")
        P,Q = parse.(Int,split(gt_s, ":")[2:3])
        KGrid(Hofstadter{P,Q}, 2, Ns, t, tp, tpp)
    else
        throw(ArgumentError("Unkown grid type: $kg"))
    end
end

function gen_shifted_ϵkGrid(kg::KGrid,shift::NTuple)  
    D = length(kg.kGrid[1])
    if D != length(shift)
        throw(ArgumentError("Grid dimension differs from shift dimension!"))
    else
        shifted_kgrid = Vector{NTuple{D,Float64}}(undef,length(kg.kGrid))
        for (i,k) in enumerate(kg.kGrid)
            shifted_kgrid[i] = k .+ shift
        end
        return kg.gen_ϵkGrid(shifted_kgrid,kg.t)
    end
end
