"""

Green's function


"""
module Greenfunc

export GreenUR
using StaticArrays, GreenFunc

struct GreenUR{T<:AbstractFloat}
    timetype::Symbol
    spacetype::Symbol

    color::Int
    timegrid::AbstractArray{T,1}
    spacegrid::AbstractArray{T,1}
    value::AbstractArray{T,3}
    function GreenUR{T}(ttype,stype,color_n,tgrid,sgrid,val)where{T<:AbstractFloat}
        @assert length(size(tgrid))==1    
        @assert length(size(sgrid))==1
        @assert length(size(val))==3
        @assert size(val)[1] == length(tgrid) && size(val)[2] == length(sgrid) && size(val)[3] == color_n
        timetype = ttype
        spacetype = stype
        color = color_n
        timegrid = tgrid
        spacegrid = sgrid
        value = val
        return new{T}(timetype,spacetype,color,timegrid,spacegrid,value)
    end 
end

Base.getindex(green::GreenUR, i ,j, k) = green.value[i,j,k]
Base.getindex(green::GreenUR, i) = green.value[i]
#Base.lastindex(green::GreenUR) = grid.size


end
#function Base.getindex
 
