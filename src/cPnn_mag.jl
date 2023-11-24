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
    kMult = ones(length(ind))
    expand_perms = map(x -> [CartesianIndex{2}(x)],ind[:])
    red_map = CartesianIndex.(ind[:])
    red_conv_map = reverse(red_map)
    return ind, red_map, red_conv_map, kMult, expand_perms, kGrid[:]
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

function gen_ϵkGrid_Hofstadter_ij(::Type{Hofstadter{P,Q}}, i::Int, j::Int, kᵢ::NTuple{2,Float64}, t::Float64, tp::Float64, tpp::Float64) where {P,Q}
    res = 0.0

    xv = kᵢ[1] + 2*π*P*(i-1)/Q
    yv = 1im*Q*kᵢ[2]
    if i==j
        res += -2*t*(cos(xv)) - 2*tpp*cos(2*xv)
    end
    if i == j+1 || i == j-1
        res += -t - 2*tp*cos(xv)
    end
    if i == j+2 || i == j-2
        res += -tpp
    end
    if i == Q && j == 1
        tp_x = kᵢ[1] - π*P/Q
        res += -t*exp(-yv) -2*tpp*cos(tp_x)*exp(-yv)
    end
    if j == Q && i == 1
        tp_x = kᵢ[1] - π*P/Q
        res += -t*exp(yv) -2*tpp*cos(tp_x)*exp(yv)
    end
    if i == Q && j == 2
        res += -tpp*exp(-yv)
    end
    if j == Q && i == 2
        res += -tp*exp(yv)
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
