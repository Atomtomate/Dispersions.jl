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
    function KGrid(GT::Type{T}, D::Int, Nk::Int, t::Float64; fftw_plan=nothing) where T<:KGridType
        sampling=[(2*π/Nk) * j - π for j in 1:Nk]
        kGrid  = collect(Base.product([sampling for Di in 1:D]...))[:]
        #TODO: = gen_kGridVecs
        #TODO: remove plan?
        fftw_plan = fftw_plan === nothing ? plan_fft(randn(Complex{Float64}, repeat([Nk], D)...), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new{GT,D}(Nk^D, Nk, kGrid, gen_ϵkGrid(SC,kGrid,t),t,Array{ComplexF64,D}(undef,repeat([Nk],D)...),fftw_plan)
    end
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #
gen_ϵkGrid(::Type{SC}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -2*t*sum(cos.(kᵢ)), kGrid))
gen_ϵkGrid(::Type{FCC}, kGrid::GridPoints, t::T) where T <: Real = collect(map(kᵢ -> -2*t*(cos(kᵢ[1])*cos(kᵢ[2])+cos(kᵢ[1])*cos(kᵢ[3])+cos(kᵢ[2])*cos(kᵢ[3])), kGrid))

gen_kGridVecs(::Type{SC}) = throw("Not implemented yet")
gen_kGridVecs(::Type{FCC}) = throw("Not implemented yet")

#TODO: match type of kGRid from KGrid{T}
ifft_post(kG::KGrid, x::Array{T,N}) where {N, T <: Number} = ShiftedArrays.circshift(x, floor.(Int, gridshape(kG) ./ 2) .+ 1)
ifft_post!(kG::KGrid, res::AbstractArray{T,D}, x::AbstractArray{T,D}) where {D, T <: Number} = ShiftedArrays.circshift!(res, x, floor.(Int, gridshape(kG) ./ 2) .+ 1)
