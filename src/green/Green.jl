"""

General Green's function. 


"""
module GreenBasic

export Green2,Green2DLR,toTau, toMatFreq, toDLR, get
using StaticArrays, Lehmann, CompositeGrids

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
- 'name': Name of green's function
- 'isFermi': Particle is fermi or boson
- 'β': Inverse temperature
- 'timeType': Whether the Green's function is in time/frequency/dlr space
- 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
- 'timeGrid': Time or Frequency grid
- 'spaceType': Whether the Green's function is in coordinate space/momentum space
- 'spaceSymmetry': Symmetry of lattice
- 'spaceGrid': Coordinate or momentum grid
- 'color': Indices of species of Green's function (such as different spin values)
- 'dlrGrid': In-built Discrete Lehmann Representation
- 'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space.
- 'dynamic': Dynamic part of Green's function
- 'instantError': Error of instantaneous part
- 'dynamicError': Error of dynamic part
"""
mutable struct Green2DLR{T<:Number,TGT,SGT,CT}
    name::Symbol
    isFermi::Bool
    β::Float64
    timeType::Symbol
    timeSymmetry::Symbol
    timeGrid::AbstractGrid
    spaceType::Symbol
    spaceSymmetry::Any
    spaceGrid::AbstractGrid
    color::CT
    dlrGrid::DLRGrid
    instant::AbstractArray{T,3}
    dynamic::AbstractArray{T,4}
    hasError::Bool
    instantError::AbstractArray{T,3}
    dynamicError::AbstractArray{T,4}

    """
         function Green2DLR{T}(timeType,spaceType,isFermi, β, Euv, rtol,spaceGrid::SGT; color::CT=nothing, timeGrid::TGT=nothing,timeSymmetry=:none, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
    
    Create two-leg Green's function on timeGrid, spaceGrid and color, with in-built DLR.
    The value and error of instant and dynamic parts are initialized with zero.
    #Arguements
    - 'name': Name of green's function. Default: :green
    - 'isFermi': Particle is fermi or boson
    - 'Euv': the UV energy scale of the spectral density
    - 'rtol': tolerance absolute error
    - 'spaceType': Whether the Green's function is in coordinate space/momentum space
    - 'spaceSymmetry': Symmetry of lattice
    - 'spaceGrid': k/x grid
    - 'β': Inverse temperature
    - 'timeType': Whether the Green's function is in time/frequency/dlr space
    - 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
    - 'timeGrid': τ/n/ω  grid. Default: DLR grid in timeType (:τ/:n/:ω) 
    - 'color': Indices of species of Green's function (such as different spin values). Default: One element array [1]
         #Required functions
         - 'getIndex(color)': Return the index of given color
    - 'dlrGrid': In-built Discrete Lehmann Representation
    - 'hasError': If Green's function is noisy
    """
    function Green2DLR{T}(isFermi, Euv, rtol, spaceType, spaceGrid::SGT, β, timeType; timeSymmetry=:none, timeGrid::TGT=nothing, color::CT=nothing, spaceSymmetry=nothing, hasError = false, name = :green)where{T<:Number,TGT,SGT,CT}
        @assert timeType == :n ||timeType == :τ ||timeType == :ω
        @assert spaceType == :k ||spaceType == :x
        @assert timeSymmetry == :ph ||timeSymmetry == :pha ||timeSymmetry == :none
        dlrGrid = DLRGrid(Euv, β, rtol, isFermi, timeSymmetry)
        if(!(TGT <: AbstractGrid))
            if(isnothing(timeGrid))
                if(timeType == :n)
                    timeGrid = dlrGrid.n
                elseif(timeType == :τ)
                    timeGrid = dlrGrid.τ
                elseif(timeType == :ω)
                    timeGrid = dlrGrid.ω
                end
                timeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(timeGrid)}(timeGrid)
            elseif(TGT <: AbstractVector{})
                timeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(timeGrid)}(timeGrid)
            else
                error("Input timeGrid has to be Vector or Composite grid")
            end
        end

        if(!(SGT <: AbstractGrid))
            if(SGT <: AbstractVector{})
                spaceGrid = CompositeGrids.SimpleG.Arbitrary{eltype(spaceGrid)}(spaceGrid)
            else
                error("Input spaceGrid has to be Vector or Composite grid")
            end
        end
     
        ct = CT
        if(isnothing(color))
            color = [1]
            ct = typeof(color)
        end
        instant = zeros(T, (length(color),length(color),spaceGrid.size))
        dynamic = zeros(T, (length(color), length(color), spaceGrid.size,timeGrid.size))
        if hasError
            instantError = zeros(T, (length(color),length(color),spaceGrid.size))
            dynamicError = zeros(T, (length(color), length(color), spaceGrid.size,timeGrid.size))
        else
            instantError = Array{T,3}(undef, 0, 0, 0)
            dynamicError = Array{T,4}(undef, 0, 0, 0, 0)
        end

        return new{T,typeof(timeGrid),typeof(spaceGrid),ct}(name, isFermi, β, timeType, timeSymmetry, timeGrid, spaceType, spaceSymmetry, spaceGrid, color, dlrGrid, instant,dynamic, hasError, instantError, dynamicError)
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
    function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
    Convert Green's function to τ space by Fourier transform.
    If green is already in τ space then it will be interpolated to the new grid.
    #Arguements
    - 'green': Original Green's function
    - 'targetGrid': Grid of outcome Green's function. Default: DLR τ grid
"""
function toTau(green::Green2DLR, targetGrid = green.dlrGrid.τ)
    if green.hasError
        #error = green.dynamicError #Need to fix Lehmann bug
        error = nothing
    else
        error = nothing
    end
    if(typeof(targetGrid) <: AbstractVector{})
        targetGrid =  CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end
    if(green.timeType == :τ)
        dynamic = tau2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis=4)
    elseif(green.timeType == :n)
        dynamic = matfreq2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis=4)
    elseif(green.timeType == :ω)
        dynamic = dlr2tau(green.dlrGrid, green.dynamic, targetGrid.grid; axis=4)
    end
    green_new = Green2DLR{eltype(dynamic)}(green.isFermi, green.dlrGrid.Euv, green.dlrGrid.rtol, green.spaceType, green.spaceGrid, green.β, :τ; timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, color = green.color, spaceSymmetry = green.spaceSymmetry)
    green_new.dynamic = dynamic
    green_new.instant = green.instant
    if green.hasError
        green_new.instantError = green.instantError
        green_new.dynamicError = green_new.dynamic
        # Need future implementation
    end
    return green_new
end

"""
    function toMatFreq(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
    Convert Green's function to ωn space by Fourier transform.
    If green is already in ωn space then it will be interpolated to the new grid.
    #Arguements
    - 'green': Original Green's function
    - 'targetGrid': Grid of outcome Green's function. Default: DLR n grid
"""
function toMatFreq(green::Green2DLR, targetGrid = green.dlrGrid.n)
    if green.hasError
        #error = green.dynamicError
        error = nothing
    else
        error = nothing
    end
    if(typeof(targetGrid) <: AbstractVector{})
        targetGrid =  CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end
    if(green.timeType == :n)
        dynamic = matfreq2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis=4)      
    elseif(green.timeType == :τ)
        dynamic = tau2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis=4)        
    elseif(green.timeType == :ω)
        dynamic = dlr2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid; axis=4)
    end
    green_new = Green2DLR{eltype(dynamic)}(green.isFermi, green.dlrGrid.Euv, green.dlrGrid.rtol, green.spaceType, green.spaceGrid, green.β, :n; timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, color = green.color, spaceSymmetry = green.spaceSymmetry)
    green_new.dynamic = dynamic
    green_new.instant = green.instant
    if green.hasError
        green_new.instantError = green.instantError
        green_new.dynamicError = green_new.dynamic
        # Need future implementation
    end
    return green_new
end

"""
    function toTau(green::Green2DLR)
    Convert Green's function to its DLR coefficents.
    Return a copy of green If it is already in dlr space.
    #Arguements
    - 'green': Original Green's function
"""
function toDLR(green::Green2DLR)
    if green.hasError
        #error = green.dynamicError
        error = nothing
    else
        error = nothing
    end
    targetGrid =  CompositeGrids.SimpleG.Arbitrary{eltype(green.dlrGrid.ω)}(green.dlrGrid.ω)
    if(green.timeType == :ω)
        green_new = deepcopy(green)
        return green_new
    elseif(green.timeType == :n)
        dynamic = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; error,axis=4)
    elseif(green.timeType == :τ)
        dynamic = tau2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; error,axis=4)
    end
    green_new = Green2DLR{eltype(dynamic)}(green.isFermi, green.dlrGrid.Euv, green.dlrGrid.rtol, green.spaceType, green.spaceGrid, green.β, :ω; timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, color = green.color, spaceSymmetry = green.spaceSymmetry)
    green_new.dynamic = dynamic
    green_new.instant = green.instant
    if green.hasError
        green_new.instantError = green.instantError
        green_new.dynamicError = green_new.dynamic
        # Need future implementation
    end
    return green_new
end


"""
    function getValue(green::Green2DLR, time, space; color=nothing, timeMethod = :default, spaceMethod=:default)

Find value of Green's function at given color, τ/ωn and k/x by interpolation.
Interpolation in τ/ωn use DLR method
#Argument
- 'green': Green's function
- 'time': Target τ/ωn point 
- 'space': Target k/x point
- 'colorvalue': Target color
- 'timeMethod':Method of interpolation in ωn/τ
- 'spacemethod': Method of interpolation in k/x 
"""
function get(green::Green2DLR, time, space; color=nothing, timeMethod=:default, spaceMethod=:default)
    @assert  green.timeType == :n ||green.timeType == :τ
    if green.hasError
        #error = green.dynamicError
        error = nothing
    else
        error = nothing
    end
    if(isnothing(color))
        cindex = 1
    else
        cindex = green.color.getIndex(color)
    end
    dynamic_c = green.dynamic[cindex, cindex, :, :]
    dynamic_x = CompositeGrids.Interp.interp1D(dynamic_c, green.spaceGrid, space; axis=1, interpstyle=spaceMethod)
    if(timeMethod == :dlr)
        if(green.timeType == :n)
            dynamic_τ = (matfreq2matfreq(green.dlrGrid, dynamic_x, [timevalue,], green.timeGrid.grid; error))[1]
        elseif(green.timeType == :τ)
            dynamic_τ = (tau2tau(green.dlrGrid, dynamic_x, [timevalue,], green.timeGrid.grid; error))[1]     
        end
    else
        dynamic_τ = CompositeGrids.Interp.interp1D(dynamic_c, green.timeGrid, time; axis=1, interpstyle=timeMethod)
    end
    return dynamic_τ
end

end

 
