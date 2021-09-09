function Base.show(io::IOT, gr::KGrid{T, D}) where {IOT <: IO, T<:KGridType, D}
    if get(io, :compact, true)
        print(io, "Reduced $T grid in $D dimensions with ", Nk(gr), " k-points.")
    else
        print(io, "Reduced $T grid in $D dimensions with ", Nk(gr), " k-points.\n
   === Public Fields ===
Nk    : total number of points
Ns    : number of points per dimension
kGrid : gridpoints
kInd  : indices
kMult : multiplicity per point
t     : hopping amplitude
ϵkGrid: dispersion")
    end
end

Base.show(io::IOT, ::MIME"text/plain", gr::KGrid{T, D}) where {IOT <: IO, T<:KGridType, D} = show(io, gr)

#= struct ReducedKIndices_cP{Ind<:Union{rInd2D, rInd3D}} <: ReducedKIndices{Array{Ind,1}} =# 
#= struct ReducedKPoints_cP{gInd<:Union{rGridP2D, rGridP3D}} <: ReducedKPoints{Array{gInd,1}} =# 
#= struct ReducedKGrid_cP{IndT <: ReducedKIndices_cP, GridT <: ReducedKPoints_cP}  <: KGrid_cP{IndT, GridT} =#
#=     Nk::Int64 =#
#=     kInd::IndT =#
#=     kMult::Array{Float64,1} =#
#=     kGrid::GridT =#
#= end =#
