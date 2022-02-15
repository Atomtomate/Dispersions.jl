import Base.collect

# ================================================================================ #
#                                Simple Cubic 2D                                    #
# ================================================================================ #

# -------------------------------------------------------------------------------- #
#                                     Types                                        #
# -------------------------------------------------------------------------------- #

abstract type p6m <: KGridType end

"""
    FullKGrid_p6m  <: FullKGrid{p6m}

Fields
-------------
- **`Nk`**     : `Int` k points in full grid.
- **`Ns`**     : `Int` k points along one dimension.
- **`kGrid`**  : `Array{Tuple{Float64, ...}}` of kGrids. Each element is a D-tuple
- **`ϵkGrid`** : `Array{Tuple{Float64, ...}}` dispersion relation.
- **`t`**      : `Float64` hopping parameter
"""
struct FullKGrid_p6m  <: FullKGrid{p6m,2}
    Nk::Int
    Ns::Int
    kGrid::GridPoints{2}
    ϵkGrid::GridDisp
    t::Float64
    fftw_plan::FFTW.cFFTWPlan
    function FullKGrid_p6m(Nk::Int, t::Float64; fftw_plan=nothing)
        ki_monkhorst = [(2*j - Nk - 1)/(2*Nk) for j in 1:Nk]
        ki_mesh = Base.product(ki_monkhorst,ki_monkhorst);
        #f_old(x) = (2π * (x[1] - x[2]/sqrt(3)), 4π*x[2]/sqrt(3))
        f(x) =  (2π * x[1], -2π * x[1]/sqrt(3) + 4π*x[2]/sqrt(3))
        kGrid = collect(map(f, ki_mesh))[:]
        disp = gen_ϵkGrid(p6m, kGrid, t)
        fftw_plan = fftw_plan === nothing ? plan_fft!(randn(Complex{Float64}, Nk, Nk), flags=FFTW.ESTIMATE, timelimit=Inf) : fftw_plan
        new(Nk^2, Nk, kGrid, disp, t, fftw_plan)
    end
end

"""
    ReducedKGrid_p6m  <: ReducedKGrid{p6m}

Reduced k grid only containing points necessary for computation of quantities involving 
this grid type and multiplicity of each stored point.

Fields
-------------
- **`Ns`**    : `Int` k points in each dimension.
- **`Nk`**    : `Int` Total number of k points.
- **`kInd`**  : `Array{Tuple{Int,...}}` indices in full grid. used for reconstruction of full grid.
- **`kMult`** : `Array{Float64,1}` multiplicity of point, used for calculations involving reduced k grids.
- **`kGrid`** : `Array{Tuple{Float64,...}}` k points of reduced grid.
"""
struct ReducedKGrid_p6m  <: ReducedKGrid{p6m, 2}
    Nk::Int
    Ns::Int
    t::Float64
    kInd::GridInd{2}
    kMult::Array{Float64,1}
    kGrid::GridPoints{2}
    ϵkGrid::GridDisp
    expand_perms::Vector{Vector{CartesianIndex{2}}}
    expand_cache::Array{Complex{Float64},2}
    fftw_plan::FFTW.cFFTWPlan
end

# -------------------------------------------------------------------------------- #
#                                   Interface                                      #
# -------------------------------------------------------------------------------- #

gridshape(kG::FullKGrid_p6m) = (kG.Ns, kG.Ns)
gridshape(kG::ReducedKGrid_p6m) = (kG.Ns, kG.Ns)

# ---------------------------- BZ to f. irr. BZ -------------------------------

"""
    reduceKGrid(kG::FullKGrid{p6m})

Returns the grid on the fully irredrucible BZ.
Filters an arbitrary grid, defined on a full kGrid, so that only 
the lower triangle remains, i.e.
for any (x_1, x_2, ...) the condition x_1 >= x_2 >= x_3 ... 
is fulfilled.
"""
function reduceKGrid(kG::FullKGrid{p6m})
    #kGrid = reshape(kG.kGrid, gridshape(kG))
    #ϵkGrid = reshape(kG.ϵkGrid, gridshape(kG))
    ind = collect(Base.product([1:kG.Ns for Di in 1:2]...))

    #nh  = floor(Int,gridshape(kG)[2]/2) + 1
    #index = ind[:,nh:end] 

    #ind_red = GridInd2D(undef, length(index))
    #grid_red = GridPoints2D(undef, length(index))
    #ϵk_red = GridDisp(undef, length(index))

    #for (i,ti) in enumerate(index)
    #    ind_red[i] = ind[ti...]
    #    grid_red[i] = kGrid[ti...]
    #    ϵk_red[i] = ϵkGrid[ti...]
    #end
	#kmult = kGrid_multiplicity(p6m, ind_red)
    #return ReducedKGrid_p6m(kG.Nk, kG.Ns, kG.t, ind_red, kmult, grid_red, ϵk_red)
    expand_cache = Array{Complex{Float64}}(undef, gridshape(kG))
    return ReducedKGrid_p6m(kG.Nk, kG.Ns, kG.t, ind[:], ones(length(ind)), kG.kGrid[:], kG.ϵkGrid[:],
                            map(x -> [CartesianIndex{2}(x)],ind[:]),expand_cache, kG.fftw_plan)
end

gen_ϵkGrid(::Type{p6m}, kGrid::GridPoints{2}, t::T1) where T1 <: Number = collect(map(kᵢ -> -2*t*(cos.(0.5*(kᵢ[1] + sqrt(3)*kᵢ[2])) + cos(0.5*(kᵢ[1] - sqrt(3)*kᵢ[2])) + cos(kᵢ[1])), kGrid))

function ifft_post(kG::ReducedKGrid_p6m, x::Array{T,2}) where T <: Number
    nh = trunc(Int, size(x,2)/2)
    x[1:nh,1:nh],x[nh+1:end,nh+1:end] = x[nh+1:end,nh+1:end],x[1:nh,1:nh]
    x[1:nh,nh+1:end],x[nh+1:end,1:nh] = x[nh+1:end,1:nh],x[1:nh,nh+1:end]
    return x
end
