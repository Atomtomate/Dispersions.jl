import Base.collect
#TODO update tests and docu for improved expandKArr


"""
    gen_SC_kGrid(Nk::Int64, D::Int64, t::Float64)

Generates a simple cubic lattice in `D` Dimensions
See also [`KGrid_SC_3D`](@ref) and [`KGrid_SC_3D`](@ref)
# Examples
```
julia> gen_SC_kGrid(2, 2)
"""
function gen_SC_kGrid(Nk::Int64, D::Int64, t::Float64)
    if D == 2 || D == 3
        return KGrid(SC, D, Nk, t)
    else
        throw("Simple Cubic only implemented for 2D and 3D")
    end
end

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

#TODO: no need for lattice types, when general basis is implemented
abstract type SC <: KGridType end
abstract type FCC <: KGridType end

# -------------------------------------------------------------------------------- #
#                                 Implementation                                   #
# -------------------------------------------------------------------------------- #

"""
    KGrid{D}  <: KGrid{SC{D}}

Fields
-------------
- **`kGrid`** : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
#TODO: Docu
"""
struct KGrid{T <: KGridType, D}
    Nk::Int
    Ns::Int
    kGrid::Array{NTuple{D,Float64},1}
    ϵkGrid::GridDisp
    t::Float64
    fft_cache::Array{ComplexF64,D}
    fftw_plan::FFTW.cFFTWPlan
    function KGrid(GT::Type{T}, D::Int, Nk::Int, t::Float64; rot_angles=nothing, fftw_plan=nothing) where T<:KGridType
        sampling=[(2*π/Nk) * j - π for j in 1:Nk]
        rot_angles = if rot_angles==nothing 
                        (D == 2) ? (0.0,) : (0.0,0.0,0.0)
                    else
                        rot_angles
                    end
        kGrid  = basis_transform(GT, Nk, collect(Base.product([sampling for Di in 1:D]...)), angles=rot_angles)[:]
        #TODO: remove plan?
        fftw_plan = fftw_plan === nothing ? plan_fft(randn(Complex{Float64}, repeat([Nk], D)...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{GT,D}(Nk^D, Nk, kGrid, gen_ϵkGrid(SC,kGrid,t),t,Array{ComplexF64,D}(undef,repeat([Nk],D)...),fftw_plan)
    end
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #
#TODO: this can be generalized for arbitrary basis vectors
gen_ϵkGrid(::Type{SC}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -2*t*sum(cos.(kᵢ)), kGrid))
gen_ϵkGrid(::Type{FCC}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -4*t*(cos(kᵢ[1]/2)*cos(kᵢ[2]/2)+cos(kᵢ[1]/2)*cos(kᵢ[3]/2)+cos(kᵢ[2]/2)*cos(kᵢ[3]/2)), kGrid))

function basis_transform(::Type{SC}, Nk::Int, kGrid::AbstractArray; angles=(0.0,0.0,0.0))
    rot = length(kGrid[1]) == 2 ?  Angle2d(angles...) : RotXYZ(angles...)
    s = 2π/Nk - π
    map(kᵢ-> Tuple( mod.(rot * collect(kᵢ) .- s, 2π) .+ s), kGrid)
end

function basis_transform(::Type{FCC}, Nk::Int, kGrid::AbstractArray; angles=(0.0,0.0,0.0))
    rot = RotXYZ(angles...)
    s = 0#2*(2π/Nk - π)
    map(kᵢ -> Tuple(rot * ([-1.0 1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 -1.0] * collect(kᵢ))), kGrid)
end

#TODO: is this the same for all k grids?
ifft_post(kG::KGrid, x::Array{T,N}) where {N, T <: Number} = reverse(x)
ifft_post!(kG::KGrid, x::AbstractArray{T,D}) where {D, T <: Number} = reverse!(x)
