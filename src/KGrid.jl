"""
    KGrid{T <: KGridType, D}

Fields
-------------
- **`Nk`**      : `Int`, Number of total k-points
- **`Ns`**      : `Int`, Number of sampling points per dimension
- **`t`**       : `Float64`, hopping parameter
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
    kGrid::GridPoints
    ϵkGrid::GridDisp
    kInd::GridInd
    kInd_conv::GridInd
    kMult::Array{Float64,1}
    expand_perms::Vector{Vector{CartesianIndex{D}}}
    cache1::Array{ComplexF64,D}
    cache2::Array{ComplexF64,D}
    fftw_plan::FFTW.cFFTWPlan
    function KGrid(GT::Type{T}, D::Int, Ns::Int, t::Float64; fftw_plan=nothing) where T<:KGridType
        sampling = gen_sampling(GT, D, Ns)
        kGrid_f = map(v -> basis_transform(GT, v), sampling)
        kInd, kInd_conv, kMult, expand_perms, kGrid = reduce_KGrid(GT, D, Ns, kGrid_f)
        ϵkGrid =  gen_ϵkGrid(GT, kGrid, t)
        gs = repeat([Ns], D)
        fftw_plan = fftw_plan === nothing ? plan_fft!(FFTW.FakeArray{ComplexF64}(gs...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{GT,D}(Ns^D, Ns, t, kGrid, ϵkGrid, kInd, kInd_conv, kMult, expand_perms,
                  Array{ComplexF64,D}(undef, gs...), Array{ComplexF64,D}(undef, gs...), fftw_plan)
    end
end

"""
    gen_kGrid(kG::String, Ns::Int)

Generates a KGrid of type and hopping strength, given in `kG` with `Ns` sampling points in the first Brillouin zone. Options are:
- '3dcP-...' : simple cubic 3D
- '2dcP-...' : simple cubic 2D
- 'cF-...'   : FCC
- 'p6m-...'  : hexagonal

# Examples
```
julia> gen_kGrid("3dcP-1.5", 10)
cP(t=1.5) grid in 3 dimensions with 1000 k-points.
```
"""
function gen_kGrid(kg::String, Ns::Int)
    findfirst("-", kg) === nothing && throw(ArgumentError("Please provide lattice type and hopping, e.g. SC3D-1.1"))
    sp = findfirst("-", kg)[1]
    data = [kg[1:(sp-1)], kg[(sp+1):end]]
    gt_s = lowercase(data[1])
    t = parse(Float64, data[2])
    if gt_s == "3dsc" || gt_s == "3dcp"
        KGrid(cP, 3, Ns, t)
    elseif gt_s == "2dsc" || gt_s == "2dcp"
        KGrid(cP, 2, Ns, t)
    elseif gt_s == "fcc" || gt_s == "cf"
        KGrid(cF, 3, Ns, t)
    elseif gt_s == "bcc" || gt_s == "ci"
        KGrid(cI, 3, Ns, t)
    elseif gt_s == "p6m"
        KGrid(p6m, 2, Ns, t)
    else
        throw(ArgumentError("Unkown grid type"))
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