# ================================================================================ #
#                            2D Cubic for Hofstadter Model                         #
# ================================================================================ #

abstract type Hofstadter{P,Q} <: KGridType end

# -------------------------------------------------------------------------------- #
#                            Grid Generator Functions                              #
# -------------------------------------------------------------------------------- #

gen_sampling(::Type{Hofstadter{P,Q}}, D::Int, Ns::Int) where {P,Q} =
    Base.product([[(2 * π / Ns) * j - π for j = 1:Ns] for Di = 1:D]...)

basis_transform(::Type{Hofstadter{P,Q}}, v::Tuple) where {P,Q} = (v[1], v[2]/Q)

function reduce_KGrid(::Type{Hofstadter{P,Q}}, D::Int, Ns::Int, kGrid::AbstractArray) where {P,Q}
    ind = collect(Base.product([1:Ns for Di in 1:2]...))
    I = LinearIndices(ind)
    index = [CartesianIndex(el) for el in ind[:]]
    kMult = ones(length(ind))
    expand_perms = map(x -> [I[CartesianIndex{2}(x)]],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    # Change from CartesianIndices to LinearIndices for performance reasons
    
    index = I[index]
    # ind_red_conv = I[ind_red_conv]
    # ind_red_crossc = I[ind_red_crossc]
    return index, red_map, red_conv_map, kMult, expand_perms, kGrid[:]
end

function gen_ϵkGrid_Hofstadter_ij_old(::Type{Hofstadter{P,Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64, tp::Float64, tpp::Float64) where {P,Q}
    res = 0.0
    if i==j
        res += -2*t*(cos(kᵢ[1] + 2*π*P*(i-1)/Q))
    end
    if i == j+1 || i == j-1
        res += -t
    end
    if i == Q && j == 1
        res += -t*exp(-1im*Q*kᵢ[2])
    end
    if j == Q && i == 1
        res += -t*exp(1im*Q*kᵢ[2])
    end
    return res
end

"""
    gen_ϵkGrid_Hofstadter_ij(::Type{Hofstadter{P,Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64, tp::Float64, tpp::Float64) where {P,Q}

Indices for entries of dispersion matrix for Hofstadter model.
"""
function gen_ϵkGrid_Hofstadter_ij(::Type{Hofstadter{P,Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64, tp::Float64, tpp::Float64) where {P,Q}
    res = 0.0

    # kᵢ[1] .=, kx, kᵢ[2] .=. ky
    # definition: epsilon[i,j], 

    # cos argument for diagonal
    xv = kᵢ[1] + 2*π*P*(i-1)/Q

    # exp arguemnt for corner
    yv = 1im*Q*kᵢ[2] 
    

    # Diagonal (contains t and t'')
    if i==j
        res += -2*t*cos(xv) 
        res += -2*tpp*cos(2*xv)
    end

    # 1. off-diagonal (contains t and t')
    if i == j+1 || i == j-1
        res += -t 
        #????
        n = max(i,j) 
        xv_tp = kᵢ[1] + 2*π*P*((n-2) + 1/2)/Q
        res += -2*tp*cos(xv_tp)
    end

    # 2. off-diagonal (contains t'')
    if i == j+2 || i == j-2
        res += -tpp
    end

    # t and t' corners
    ## (bottom left)
    if i == Q && j == 1
        tp_x = kᵢ[1] - π*P/Q
        res += -t*exp(-yv)
        res += -2*tp*cos(tp_x)*exp(-yv)
    end

    ## (top right)
    if i == 1 && j == Q
        tp_x = kᵢ[1] - π*P/Q
        res += -t*exp(yv)
        res += -2*tp*cos(tp_x)*exp(yv)
    end

    # t'' corners
    ## (bottom left)
    if i == Q-1 && j == 1
        res += -tpp*exp(-yv)
    end
    if i == Q   && j == 2
        res += -tpp*exp(-yv)
    end

    ## (top right)
    if i == 1   && j == Q-1
        res += -tpp*exp(yv)
    end
    if i == 2   && j == Q
        res += -tpp*exp(yv)
    end
    return res
end

function gen_ϵkGrid(::Type{Hofstadter{P,Q}}, kGrid::GridPoints, t::Float64, tp::Float64, tpp::Float64) where {P,Q}
    if !(P==1 && Q==1) && (gcd(P,Q) != 1 || Q < 3)
        throw(ArgumentError("gcd(P,Q) = $(gcd(P,Q)) != 1 or Q = $(Q) < 3!"))
    end
    res = Array{ComplexF64,3}(undef, Q, Q, length(kGrid))
    for (l, kᵢ) in enumerate(kGrid)
        for j in 1:Q
            for i in 1:Q
                res[i,j,l] = gen_ϵkGrid_Hofstadter_ij(Hofstadter{P,Q}, i, j, kᵢ, t, tp, tpp)
            end
        end
    end
    res
end

# -------------------------------------------------------------------------------- #
#                             Custom Helper Functions                              #
# -------------------------------------------------------------------------------- #

conv_sample_post(kG::KGrid{Hofstadter{P,Q}, 2}, x) where {P,Q} = x
