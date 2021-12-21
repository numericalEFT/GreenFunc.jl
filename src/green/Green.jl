"""

General Green's function. 


"""
module GreenBasic

export Green2,Green2DLR,TimeFourier
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
    instant::AbstractArray{T,3}
    dynamic::AbstractArray{T,4}


    """
         function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}

    create two-leg Green's function on tgrid, sgrid and color_n.
    The value of Green's function is initialized with zero.
    """
    function Green2{T}(timeType,spaceType,isFermi, β,timeGrid::TGT, spaceGrid::SGT; color::CT=nothing,timeSymmetry=:none, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
        ct = CT
        if(color == nothing)
            color = [1]
            ct = typeof(color)
        end
        static_val = zeros(T, (length(color),length(color),length(spaceGrid)))
        dynamic_val = zeros(T, (length(color), length(color), length(spaceGrid),length(timeGrid)))
        instant = static_val
        dynamic = dynamic_val
        return new{T,TGT,SGT,ct}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry, β,color,timeGrid,spaceGrid,instant,dynamic, error)
    end
end


"""
Green's function with two external legs that has in-built Discrete Lehmann Representation.

#Members:
- 'timeType': Whether the Green's function is in time/frequency/dlr space
- 'spaceType': Whether the Green's function is in coordinate space/momentum space
- 'isFermi': Particle is fermi or boson
- 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
- 'spaceSymmetry': Symmetry of lattice
- 'β': Inverse temperature
- 'Euv': the UV energy scale of the spectral density
- 'rtol': tolerance absolute error
- 'color': Indices of species of Green's function (such as different spin values)
- 'timeGrid': Time or Frequency grid
- 'spaceGrid': Coordinate or momentum grid
- 'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space.
- 'dynamic': Dynamic part of Green's function
- 'error': The error of noisy Green's function
- 'dlrGrid': In-built Discrete Lehmann Representation
"""

mutable struct Green2DLR{T<:Number,TGT,SGT,CT}
    timeType::Symbol
    spaceType::Symbol
    isFermi::Bool
    timeSymmetry::Symbol
    spaceSymmetry::Any
    β::Float64
    Euv::Float64
    rtol::Float64
    color::CT
    timeGrid::TGT
    spaceGrid::SGT
    instant::AbstractArray{T,3}
    dynamic::AbstractArray{T,4}
    error::Any
    dlrGrid::DLRGrid

    """
         function Green2DLR{T}(timeType,spaceType,isFermi, β, Euv, rtol,spaceGrid::SGT; color::CT=nothing, timeGrid::TGT=nothing,timeSymmetry=:none, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
    
    Create two-leg Green's function on timeGrid, spaceGrid and color, with in-built DLR.
    The value of instant and dynamic parts are initialized with zero.
    #Arguements
    - 'timeType': Whether the Green's function is in time/frequency/dlr space
    - 'spaceType': Whether the Green's function is in coordinate space/momentum space
    - 'isFermi': Particle is fermi or boson
    - 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
    - 'spaceSymmetry': Symmetry of lattice
    - 'β': Inverse temperature
    - 'Euv': the UV energy scale of the spectral density
    - 'rtol': tolerance absolute error
    - 'color': Indices of species of Green's function (such as different spin values). Default: One element array [1]
    - 'timeGrid': Time or Frequency grid. Default: DLR grid in timeType (frequency, tau or dlr) 
    - 'spaceGrid': Coordinate or momentum grid
    - 'error': The error of noisy Green's function
    - 'dlrGrid': In-built Discrete Lehmann Representation
    """
    function Green2DLR{T}(timeType,spaceType,isFermi, β, Euv, rtol,spaceGrid::SGT; color::CT=nothing, timeGrid::TGT=nothing,timeSymmetry=:none, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
        @assert timeType == :freq ||timeType == :tau ||timeType == :dlr
        @assert spaceType == :mom ||timeType == :spa
        @assert timeSymmetry == :ph ||timeSymmetry == :pha ||timeSymmetry == :none
        dlrGrid = DLRGrid(Euv, β, rtol, isFermi, timeSymmetry)
        ct = CT
        tgt = TGT
        if(timeGrid == nothing)
            if(timeType == :freq)
                timeGrid = dlrGrid.n
            elseif(timeType == :tau)
                timeGrid = dlrGrid.τ
            elseif(timeType == :dlr)
                timeGrid = dlrGrid.ω
            end
            tgt = typeof(timeGrid)
        end
        if(color == nothing)
            color = [1]
            ct = typeof(color)
        end
        static_val = zeros(T, (length(color),length(color),length(spaceGrid)))
        dynamic_val = zeros(T, (length(color), length(color), length(spaceGrid),length(timeGrid)))
        instant = static_val
        dynamic = dynamic_val
        return new{T,tgt,SGT,ct}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry, β, Euv, rtol, color,timeGrid,spaceGrid,instant,dynamic, error,dlrGrid)
    end

    # function Green2DLR{T}(timeType,spaceType, dlrGrid::DLRGrid, color::CT,timeGrid::TGT,spaceGrid::SGT; instant=nothing, dynamic = nothing, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
    #     Euv = dlrGrid.Euv
    #     β = dlrGrid.β
    #     timeSymmetry = dlrGrid.symmetry
    #     isFermi = dlrGrid.isFermi
    #     rtol = dlrGrid.rtol
    #     if(instant == nothing)
    #         instant_val = zeros(T, (length(color),length(spaceGrid)))
    #         instant = instant_val
    #     end
    #     if(dynamic == nothing)
    #         dynamic_val = zeros(T, (length(color), length(color), length(spaceGrid),length(timeGrid)))
    #         dynamic = dynamic_val            
    #     end
    #     return new{T,TGT,SGT,CT}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry,β,Euv,rtol,color,timeGrid,spaceGrid,instant,dynamic, error,dlrGrid)
    # end
end

"""
    function TimeFourier(green::Green2DLR, aimType::Symbol; aimGrid = nothing)
    Convert Green's function to different time representation (frequency, tau , or dlr).
    #Arguements
    - 'green': Original Green's function
    - 'aimType': Time representation of outcome Green's function
    - 'aimGrid': Grid of outcome Green's function. Default: DLR grid of aimType
"""


function TimeFourier(green::Green2DLR, aimType::Symbol; aimGrid = nothing)
    if(aimGrid == nothing)
        if(aimType == :freq)
            aimGrid = green.dlrGrid.n
        elseif(aimType == :tau)
            aimGrid = green.dlrGrid.τ
        elseif(aimType == :dlr)
            aimGrid = green.dlrGrid.ω
        end
    end
    if(green.timeType==:freq)
        if(aimType == :tau)
            data_mid = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid; green.error,axis=4)
            data = dlr2tau(green.dlrGrid, data_mid, aimGrid; axis=4)
        elseif(aimType == :dlr)
            data = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid; green.error,axis=4)
        end
    elseif(green.timeType == :tau)
        if(aimType == :freq)
            data_mid = tau2dlr(green.dlrGrid, green.dynamic, green.timeGrid; green.error,axis=4)
            data = dlr2matfreq(green.dlrGrid, data_mid, aimGrid; axis=4)
        elseif(aimType == :dlr)
            data = tau2dlr(green.dlrGrid, green.dynamic, green.timeGrid; green.error,axis=4)
        end
    elseif(green.timeType == :dlr)
        if(aimType == :freq)
            data = dlr2matfreq(green.dlrGrid, green.dynamic, aimGrid; axis=4)
        elseif(aimType == :tau)
            data = dlr2tau(green.dlrGrid, green.dynamic, aimGrid; axis=4)
        end
    end

    green_new = Green2DLR{typeof(data[1,1,1,1])}(aimType,green.spaceType, green.isFermi, green.β, green.Euv, green.rtol, green.spaceGrid; color = green.color, timeGrid = aimGrid, timeSymmetry = green.timeSymmetry, spaceSymmetry = green.spaceSymmetry, error = green.error)
    green_new.dynamic = data
    green_new.instant = green.instant
    return green_new
end

end

 
