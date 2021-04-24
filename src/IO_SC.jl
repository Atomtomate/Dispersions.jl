#Base.show(io::IO, fi::Union{FullKIndices_cP,FullKPoints_cP}) = 
#        print(io, typeof(fi) <: FullKIndices ? "Index" : "GridPoints", "[",minimum(fi.grid),
#              ":",maximum(fi.grid),"]: length:", length(fi.grid))

#function Base.show(io::IO, ::MIME"text/plain", fi::Union{FullKIndices_cP{T},FullKPoints_cP{T}}) where T
#    typestr = typeof(fi) <: FullKIndices ? "Index" : "GridPoints"
#    if get(io, :compact, true)
#        print(io, typestr, "[",minimum(fi.grid),":",maximum(fi.grid),"]: length:", length(fi.grid))
#    else
#        print(io, typestr, " for full k simple cubic grid in ", T == Ind2D ? "2" : "3",
#              " Dimensions. From ", minimum(fi.grid), " to ", maximum(fi.grid), " with ", 
#              length(fi.grid), " points.")
#    end
#end

#Base.show(io::IO, gr::FullKGrid_cP{T1,T2}) where {T1,T2} = print(io, "FullKGrid_cP[", Nk(gr), "] for ", 
#                                                                 T1 <: FullKIndices_cP{Ind2D} ? "2" : "3" , "D")

#function Base.show(io::IO, ::MIME"text/plain", gr::FullKGrid_cP{T1,T2}) where {T1,T2}
#    if get(io, :compact, true)
#        print(io, "FullKGrid_cP[", Nk(gr), "] for ", T1 <: FullKIndices_cP{Ind2D} ? "2" : "3" , "D")
#    else
#        print(io, "FullKGrid for simple cubic lattice with size ", Nk(gr), " for ", 
#              T1 <: FullKIndices_cP{Ind2D} ? "2" : "3" , "D")
#    end
#end

#= struct ReducedKIndices_cP{Ind<:Union{rInd2D, rInd3D}} <: ReducedKIndices{Array{Ind,1}} =# 
#= struct ReducedKPoints_cP{gInd<:Union{rGridP2D, rGridP3D}} <: ReducedKPoints{Array{gInd,1}} =# 
#= struct ReducedKGrid_cP{IndT <: ReducedKIndices_cP, GridT <: ReducedKPoints_cP}  <: KGrid_cP{IndT, GridT} =#
#=     Nk::Int64 =#
#=     kInd::IndT =#
#=     kMult::Array{Float64,1} =#
#=     kGrid::GridT =#
#= end =#
