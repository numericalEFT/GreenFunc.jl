module MeshGrids

using ..GreenFunc
using ..Lehmann
using ..CompositeGrids

locate(m::AbstractGrid, pos) = CompositeGrids.Interp.locate(m, pos)
volume(m::AbstractGrid, index) = CompositeGrids.Interp.volume(m, index)
volume(m::AbstractGrid) = CompositeGrids.Interp.volume(m)

export locate, volume

abstract type TemporalGrid{T} <: AbstractGrid{T} end

const FERMION = true
const BOSON = false
export FERMION, BOSON

export TemporalGrid
include("common.jl")

include("imtime.jl")
export ImTime

include("imfreq.jl")
export ImFreq

include("dlrfreq.jl")
export DLRFreq

include("MeshProduct.jl")
export MeshProduct

#TODO: more functions from CompositeGrids


end
