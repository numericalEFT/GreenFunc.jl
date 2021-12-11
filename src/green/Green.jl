"""

General Green's function. 


"""
module GreenBasic

export GreenTwo
using StaticArrays, GreenFunc

"""
Green's function with two external legs. The structure saves a function G( τ, q, σ)
and corresponding grids of τ, q and σ.

#Members:
- 'timeType': Denote whether the Green's function is in time space or frequency space.
- 'spaceType': Denote whether the Green's function is in coordinate space or momentum space.
- 'color': Indices that denote different species of Green's function (such as different spin values)
- 'timeGrid': Time or Frequency grid
- 'spaceGrid': Coordinate or momentum grid
- 'value': Data of the Green's fucntion G( τ, q, σ)
"""
mutable struct GreenTwo{T<:AbstractFloat,TGT,SGT,CT}
    timeType::Symbol
    spaceType::Symbol
    color::CT
    timeGrid::TGT
    spaceGrid::SGT
    value::AbstractArray{T,3}



    """
         function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}

    create two-leg Green's function on tgrid, sgrid and color_n.
    The value of Green's function is initialized with zero.
    """

    function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}
        val = zeros(T, (length(tgrid), length(sgrid),length(color_n)))
        timeType = ttype
        spaceType = stype
        color = color_n
        timeGrid = tgrid
        spaceGrid = sgrid
        value = val
        return new{T,TGT,SGT,CT}(timeType,spaceType,color,timeGrid,spaceGrid,value)
    end 
end


Base.getindex(green::GreenTwo, i ,j, k) = green.value[i,j,k]
Base.getindex(green::GreenTwo, i) = green.value[i]
#Base.lastindex(green::GreenUR) = grid.size


end
#function Base.getindex
 
