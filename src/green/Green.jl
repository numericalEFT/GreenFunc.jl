"""

General Green's function. 


"""
module GreenBasic

export Green2
using StaticArrays, GreenFunc, Lehmann

"""
Green's function with two external legs. The structure saves a function G( τ, q, σ)
and corresponding grids of τ, q and σ.

#Members:
- 'timeType': Whether the Green's function is in time/frequency/dlr space
- 'spaceType': Whether the Green's function is in coordinate space or momentum space
- 'color': Indices of species of Green's function (such as different spin values)
- 'isFermi': Particle is fermi or boson
- 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
- 'spaceSymmetry': Symmetry of lattice
- 'timeGrid': Time or Frequency grid
- 'β':Inverse temperature
- 'spaceGrid': Coordinate or momentum grid
- 'instant': Instantaneous part of Green's function that is δ(τ) in τ space.
- 'dynamic': Dynamic part of Green's function
"""

mutable struct Green2{T<:AbstractFloat,TGT,SGT,CT}
    timeType::Symbol
    spaceType::Symbol
    isFermi::Bool
    timeSymmetry::Symbol
    spaceSymmetry::Any
    β::Float64
    color::CT
    timeGrid::TGT
    spaceGrid::SGT
    instant::AbstractArray{T,2}
    dynamic::AbstractArray{T,4}



    """
         function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}

    create two-leg Green's function on tgrid, sgrid and color_n.
    The value of Green's function is initialized with zero.
    """
    function Green2{T}(ttype,stype,isfermi,tsym, ssym, beta, color_n::CT,tgrid::TGT,sgrid::SGT)where{T<:AbstractFloat,TGT,SGT,CT}
        static_val = zeros(T, (length(color_n),length(sgrid)))        
        dynamic_val = zeros(T, (length(color_n), length(color_n), length(sgrid),length(tgrid)))
        timeType = ttype
        spaceType = stype
        isFermi = isfermi
        timeSymmetry = tsym
        spaceSymmetry = ssym 
        β = beta
        color = color_n
        timeGrid = tgrid
        spaceGrid = sgrid
        instant = static_val
        dynamic = dynamic_val
        return new{T,TGT,SGT,CT}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry,β,color,timeGrid,spaceGrid,instant,dynamic)
    end 
end

"""
Green's function with two external legs. The structure saves a function G( τ, q, σ)
and corresponding grids of τ, q and σ.

#Members:
- 'timeType': Whether the Green's function is in time/frequency/dlr space
- 'spaceType': Whether the Green's function is in coordinate space or momentum space
- 'color': Indices of species of Green's function (such as different spin values)
- 'isFermi': Particle is fermi or boson
- 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
- 'spaceSymmetry': Symmetry of lattice
- 'timeGrid': Time or Frequency grid
- 'β':Inverse temperature
- 'spaceGrid': Coordinate or momentum grid
- 'instant': Instantaneous part of Green's function that is δ(τ) in τ space.
- 'dynamic': Dynamic part of Green's function
"""

mutable struct Green2DLR{T<:AbstractFloat,TGT,SGT,CT}
    timeType::Symbol
    spaceType::Symbol
    isFermi::Bool
    timeSymmetry::Symbol
    spaceSymmetry::Any
    β::Float64
    rtol::Float64
    color::CT
    timeGrid::TGT
    spaceGrid::SGT
    instant::AbstractArray{T,2}
    dynamic::AbstractArray{T,4}
    dlrGrid::DLRGrid


    """
         function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}

    create two-leg Green's function on tgrid, sgrid and color_n.
    The value of Green's function is initialized with zero.
    """
    function Green2DLR{T}(ttype,stype,isfermi,tsym, ssym, beta, color_n::CT,tgrid::TGT,sgrid::SGT)where{T<:AbstractFloat,TGT,SGT,CT}
        static_val = zeros(T, (length(color_n),length(sgrid)))        
        dynamic_val = zeros(T, (length(color_n), length(color_n), length(sgrid),length(tgrid)))
        timeType = ttype
        spaceType = stype
        isFermi = isfermi
        timeSymmetry = tsym
        spaceSymmetry = ssym 
        β = beta
        color = color_n
        timeGrid = tgrid
        spaceGrid = sgrid
        instant = static_val
        dynamic = dynamic_val
        return new{T,TGT,SGT,CT}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry,β,color,timeGrid,spaceGrid,instant,dynamic)
    end 
end







end
#function Base.getindex
 
