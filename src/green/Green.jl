"""
General Green's function. 
"""


# """
# Green's function with two external legs. The structure saves a function G( τ, q, σ)
# and corresponding grids of τ, q and σ.

# #Members:
# - 'timeType': Whether the Green's function is in time/frequency/dlr space
# - 'spaceType': Whether the Green's function is in coordinate space or momentum space
# - 'color': Indices of species of Green's function (such as different spin values)
# - 'isFermi': Particle is fermi or boson
# - 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
# - 'spaceSymmetry': Symmetry of lattice
# - 'timeGrid': Time or Frequency grid
# - 'β':Inverse temperature
# - 'spaceGrid': Coordinate or momentum grid
# - 'instant': Instantaneous part of Green's function that is δ(τ) in τ space.
# - 'dynamic': Dynamic part of Green's function
# """
# mutable struct Green2{T<:AbstractFloat,TGT,SGT,CT}
#     timeType::Symbol
#     spaceType::Symbol
#     isFermi::Bool
#     timeSymmetry::Symbol
#     spaceSymmetry::Any
#     β::Float64
#     color::CT
#     timeGrid::TGT
#     spaceGrid::SGT
#     instant::AbstractArray{T,3}
#     dynamic::AbstractArray{T,4}


#     """
#          function GreenTwo{T,TGT,SGT,CT}(ttype,stype,color_n,tgrid,sgrid)where{T<:AbstractFloat,TGT,SGT,CT}

#     create two-leg Green's function on tgrid, sgrid and color_n.
#     The value of Green's function is initialized with zero.
#     """
#     function Green2{T}(timeType,spaceType,isFermi, β,timeGrid::TGT, spaceGrid::SGT; color::CT=nothing,timeSymmetry=:none, spaceSymmetry=nothing,error=nothing)where{T<:Number,TGT,SGT,CT}
#         ct = CT
#         if(color == nothing)
#             color = [1]
#             ct = typeof(color)
#         end
#         static_val = zeros(T, (length(color),length(color),length(spaceGrid)))
#         dynamic_val = zeros(T, (length(color), length(color), length(spaceGrid),length(timeGrid)))
#         instant = static_val
#         dynamic = dynamic_val
#         return new{T,TGT,SGT,ct}(timeType,spaceType,isFermi,timeSymmetry,spaceSymmetry, β,color,timeGrid,spaceGrid,instant,dynamic, error)
#     end
# end


"""
Green's function with two external legs that has in-built Discrete Lehmann Representation.

#Members:
- 'name': Name of green's function
- 'isFermi': Particle is fermi or boson
- 'β': Inverse temperature
- 'color': Number of different species of Green's function (such as different spin values)
- 'dlrGrid': In-built Discrete Lehmann Representation
- 'timeType': Whether the Green's function is in time/frequency/dlr space
- 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
- 'timeGrid': Time or Frequency grid
- 'spaceType': Whether the Green's function is in coordinate space/momentum space
- 'spaceGrid': Coordinate or momentum grid
- 'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space.
- 'dynamic': Dynamic part of Green's function
- 'instantError': Error of instantaneous part
- 'dynamicError': Error of dynamic part
"""
mutable struct Green2DLR{T<:Number,Type<:TimeDomain,TGT,SGT}
    name::Symbol
    isFermi::Bool
    β::Float64
    color::Int
    dlrGrid::DLRGrid

    #########    Mesh   ##############

    timeType::DataType
    timeSymmetry::Symbol
    timeGrid::TGT

    spaceType::Symbol
    spaceGrid::SGT

    ###########     data   ###########
    instant::Array{T,3}
    dynamic::Array{T,4}

    ####### statistical error handling #####
    instantError::Array{T,3}
    dynamicError::Array{T,4}

    """
        function Green2DLR{T,TimeType}(name::Symbol, β, isFermi::Bool, Euv, spaceGrid, color::Int = 1;
            timeSymmetry::Symbol = :none, rtol = 1e-8, kwargs...
        ) where {T<:Number,TimeType<:TimeDomain}
    Create two-leg Green's function on timeGrid, spaceGrid and color, with in-built DLR.
    #Arguements
    - 'name': Name of green's function. 
    - 'β': Inverse temperature
    - 'isFermi': Particle is fermi or boson
    - 'Euv': the UV energy scale of the spectral density
    - 'spaceGrid': k/x grid
    - 'color': Number of different species of Green's function (such as different spin values). Default: 1 color
    - 'timeSymmetry': Whether the Green's function has particle-hole symmetry, anti-particle-hole symmetry or none of them
    - 'rtol': tolerance absolute error
    #Optional Arguments
    - 'timeGrid': τ/n/ω  grid. Default: DLR grid in timeType (:τ/:n/:ω) 
    - 'instant': Instantaneous part of Green's function that is proportional to δ(τ) in τ space. Default: 0 everywhere
    - 'dynamic': Dynamic part of Green's function. Default: 0 everywhere
    - 'instantError': Error of instantaneous part. Default: 0 everywhere
    - 'dynamicError': Error of dynamic part. Default: 0 everywhere
    """
    function Green2DLR{T,TimeType}(name::Symbol, β, isFermi::Bool, Euv, spaceGrid, color::Int = 1;
        timeSymmetry::Symbol = :none, rtol = 1e-8, kwargs...
    ) where {T<:Number,TimeType<:TimeDomain}
        # @assert spaceType == :k || spaceType == :x
        @assert timeSymmetry == :ph || timeSymmetry == :pha || timeSymmetry == :none

        spaceType = :k #TODO: replace it with spaceGrid.type after we have a spaceGrid package

        dlrGrid = DLRGrid(Euv, β, rtol, isFermi, timeSymmetry)

        println(keys(kwargs))

        if :timeGrid in keys(kwargs)
            timeGrid = kwargs[:timeGrid]
        else
            if TimeType == ImFreq
                timeGrid = dlrGrid.n
            elseif TimeType == ImTime
                timeGrid = dlrGrid.τ
            elseif TimeType == DLRFreq
                timeGrid = dlrGrid.ω
            else
                error("$TimeType is not supported!")
            end
        end
        #if timeGrid is nothing, then set it to be the DLR grid, which is a vector of Integer or Float64
        if timeGrid isa AbstractVector
            timeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(timeGrid)}(timeGrid)
        end
        @assert timeGrid isa AbstractGrid "Input timeGrid has to be Vector or Composite grid"
        if TimeType == DLRFreq
            #@assert length(timeGrid.grid) == dlrGrid.size "The size of the DLR grid should match the DLR rank = $(dlrGrid.size)."
        elseif TimeType == ImFreq
            @assert eltype(timeGrid.grid) <: Int "Matsubara frequency grid is expected to be integers!"
        end

        # println(TGT)

        if spaceGrid isa AbstractVector
            spaceGrid = CompositeGrids.SimpleG.Arbitrary{eltype(spaceGrid)}(spaceGrid)
        end
        @assert spaceGrid isa AbstractGrid "Input spaceGrid has to be Vector or Composite grid"

        instant = Array{T,3}(undef, 0, 0, 0)
        dynamic = Array{T,4}(undef, 0, 0, 0, 0)
        instantError = Array{T,3}(undef, 0, 0, 0)
        dynamicError = Array{T,4}(undef, 0, 0, 0, 0)

        gnew = new{T,TimeType,typeof(timeGrid),typeof(spaceGrid)}(
            name, isFermi, β, color, dlrGrid,
            TimeType, timeSymmetry, timeGrid,
            spaceType, spaceGrid,
            instant, dynamic,
            instantError, dynamicError)

        return set!(gnew; kwargs...)
    end
end

function Base.size(green::Green2DLR)
    return (green.color, green.color, size(green.spaceGrid), size(green.timeGrid))
end

function set!(green::Green2DLR; kwargs...)
    dynamicSize = size(green)
    instantSize = size(green)[1:3]
    if :dynamic in keys(kwargs) && isempty(kwargs[:dynamic]) == false
        green.dynamic = reshape(kwargs[:dynamic], Tuple(dynamicSize))
    end
    if :instant in keys(kwargs) && isempty(kwargs[:instant]) == false
        green.instant = reshape(kwargs[:instant], Tuple(instantSize))
    end
    if :dynamicError in keys(kwargs) && isempty(kwargs[:dynamicError]) == false
        green.dynamicError = reshape(kwargs[:dynamicError], Tuple(dynamicSize))
    end
    if :instantError in keys(kwargs) && isempty(kwargs[:instantError]) == false
        green.instantError = reshape(kwargs[:instantError], Tuple(instantSize))
    end
    return green
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

    if targetGrid isa AbstractVector
        targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end

    # do nothing if the domain and the grid remain the same
    if green.timeType == ImTime && length(green.timeGrid.grid) ≈ length(targetGrid.grid) && green.timeGrid.grid ≈ targetGrid.grid
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImTime)
        dynamic = tau2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == ImFreq)
        dynamic = matfreq2tau(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == DLRFreq)
        dynamic = dlr2tau(green.dlrGrid, green.dynamic, targetGrid.grid; axis = 4)
    end

    return Green2DLR{eltype(dynamic),ImTime}(
        green.name, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)

end

"""
    function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
    Convert Green's function to τ space by Fourier transform.
    If green is already in τ space then it will be interpolated to the new grid.
    #Arguements
    - 'green': Original Green's function
    - 'targetGrid': Grid of outcome Green's function. Default: DLR τ grid
"""
function toMatFreq(green::Green2DLR, targetGrid = green.dlrGrid.n)

    if targetGrid isa AbstractVector
        targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
    end

    # do nothing if the domain and the grid remain the same
    if green.timeType == ImFreq && length(green.timeGrid.grid) ≈ length(targetGrid.grid) && green.timeGrid.grid ≈ targetGrid.grid
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImFreq)
        dynamic = matfreq2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == ImTime)
        dynamic = tau2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; axis = 4)
    elseif (green.timeType == DLRFreq)
        dynamic = dlr2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid; axis = 4)
    end

    return Green2DLR{eltype(dynamic),ImFreq}(
        green.name, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)

end

"""
    function toTau(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
    Convert Green's function to τ space by Fourier transform.
    If green is already in τ space then it will be interpolated to the new grid.
    #Arguements
    - 'green': Original Green's function
    - 'targetGrid': Grid of outcome Green's function. Default: DLR τ grid
"""
function toDLR(green::Green2DLR)

    
    targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(green.dlrGrid.ω)}(green.dlrGrid.ω)


    # do nothing if the domain and the grid remain the same
    if green.timeType == DLRFreq
        return green
    end
    if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
        return green
    end


    if (green.timeType == ImTime)
        dynamic = tau2dlr(green.dlrGrid, green.dynamic,  green.timeGrid.grid ; axis = 4)
    elseif (green.timeType == ImFreq)
        dynamic = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid ; axis = 4)
    end

    return Green2DLR{eltype(dynamic),DLRFreq}(
        green.name, green.β, green.isFermi, green.dlrGrid.Euv, green.spaceGrid, green.color;
        timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, rtol = green.dlrGrid.rtol,
        dynamic = dynamic, instant = green.instant)

end


# """
#     function toMatFreq(green::Green2DLR, targetGrid =  green.dlrGrid.τ)
#     Convert Green's function to ωn space by Fourier transform.
#     If green is already in ωn space then it will be interpolated to the new grid.
#     #Arguements
#     - 'green': Original Green's function
#     - 'targetGrid': Grid of outcome Green's function. Default: DLR n grid
# """
# function toMatFreq(green::Green2DLR, targetGrid = green.dlrGrid.n)
#     if green.hasError
#         #error = green.dynamicError
#         error = nothing
#     else
#         error = nothing
#     end
#     if (typeof(targetGrid) <: AbstractVector{})
#         targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
#     end
#     if (green.timeType == :n)
#         dynamic = matfreq2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis = 4)
#     elseif (green.timeType == :τ)
#         dynamic = tau2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid, green.timeGrid.grid; error, axis = 4)
#     elseif (green.timeType == :ω)
#         dynamic = dlr2matfreq(green.dlrGrid, green.dynamic, targetGrid.grid; axis = 4)
#     end
#     green_new = Green2DLR{eltype(dynamic)}(green.isFermi, green.dlrGrid.Euv, green.dlrGrid.rtol, green.spaceType, green.spaceGrid, green.β, :n; timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, color = green.color, spaceSymmetry = green.spaceSymmetry)
#     green_new.dynamic = dynamic
#     green_new.instant = green.instant
#     if green.hasError
#         green_new.instantError = green.instantError
#         green_new.dynamicError = green_new.dynamic
#         # Need future implementation
#     end
#     return green_new
# end

# """
#     function toTau(green::Green2DLR)
#     Convert Green's function to its DLR coefficents.
#     Return a copy of green If it is already in dlr space.
#     #Arguements
#     - 'green': Original Green's function
# """
# function toDLR(green::Green2DLR)
#     if green.hasError
#         #error = green.dynamicError
#         error = nothing
#     else
#         error = nothing
#     end
#     targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(green.dlrGrid.ω)}(green.dlrGrid.ω)
#     if (green.timeType == :ω)
#         green_new = deepcopy(green)
#         return green_new
#     elseif (green.timeType == :n)
#         dynamic = matfreq2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; error, axis = 4)
#     elseif (green.timeType == :τ)
#         dynamic = tau2dlr(green.dlrGrid, green.dynamic, green.timeGrid.grid; error, axis = 4)
#     end
#     green_new = Green2DLR{eltype(dynamic)}(green.isFermi, green.dlrGrid.Euv, green.dlrGrid.rtol, green.spaceType, green.spaceGrid, green.β, :ω; timeSymmetry = green.timeSymmetry, timeGrid = targetGrid, color = green.color, spaceSymmetry = green.spaceSymmetry)
#     green_new.dynamic = dynamic
#     green_new.instant = green.instant
#     if green.hasError
#         green_new.instantError = green.instantError
#         green_new.dynamicError = green_new.dynamic
#         # Need future implementation
#     end
#     return green_new
# end


"""
    function get(green::Green2DLR, time, space; color=nothing, timeMethod = :default, spaceMethod=:default)

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
function getDynamic(green::Green2DLR, time, space; colorIndex = 1, timeMethod = CompositeGrids.Interp.LinearInterp(), spaceMethod = CompositeGrids.Interp.LinearInterp())
    #@assert green.timeType == :n || green.timeType == :τ
    if !(typeof(spaceMethod)<:CompositeGrids.Interp.InterpStyle)
        error("Space Interpolation method not implemented!")
    end
    if isempty(green.dynamic)
        error("Dynamic Green's function can not be empty!")
    else
        dynamic_τx = green.dynamic[colorIndex, colorIndex, :, :]
        spaceNeighbor = CompositeGrids.Interp.findneighbor(spaceMethod, green.spaceGrid, space)
        dynamic_x_slice = CompositeGrids.Interp.dataslice(dynamic_τx, 1, spaceNeighbor.index)
        if timeMethod == DLRInterp()
            dynamic_τ = CompositeGrids.Interp.interpsliced(spaceNeighbor, dynamic_x_slice, axis = 1)
            if green.timeType == ImFreq
                dynamic = (matfreq2matfreq(green.dlrGrid, dynamic_τ, [time,], green.timeGrid.grid))[1]
            elseif green.timeType == ImTime
                dynamic = (tau2tau(green.dlrGrid, dynamic_τ, [time,], green.timeGrid.grid))[1]
            end
        elseif typeof(timeMethod)<:CompositeGrids.Interp.InterpStyle
            timeNeighbor = CompositeGrids.Interp.findneighbor(timeMethod, green.timeGrid, time)
            dynamic_τx_slice =  CompositeGrids.Interp.dataslice(dynamic_x_slice, 2, timeNeighbor.index)
            dynamic_τ = CompositeGrids.Interp.interpsliced(spaceNeighbor, dynamic_τx_slice, axis = 1)
            dynamic = CompositeGrids.Interp.interpsliced(timeNeighbor, dynamic_τ, axis = 1)
        else
            error("Time interpolation method not implemented!")
        end
    end


    return dynamic
end

function getInstant(green::Green2DLR, space; colorIndex = 1, spaceMethod = CompositeGrids.Interp.LinearInterp())
    if isempty(green.instant)
        error("Instant Green's function can not be empty!")
    elseif typeof(spaceMethod)<:CompositeGrids.Interp.InterpStyle
        instant_x = green.instant[colorIndex, colorIndex, :]
        spaceNeighbor = CompositeGrids.Interp.findneighbor(spaceMethod, green.spaceGrid, space)
        instant_x_slice = CompositeGrids.Interp.dataslice(instant_x, spaceNeighbor.index)
        instant = CompositeGrids.Interp.interpsliced(spaceNeighbor, instant_x_slice)
    else
        error("Space Interpolation method not implemented!")
    end
    return instant
end
