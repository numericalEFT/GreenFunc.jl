module MeshGrids

using ..GreenFunc
using ..Lehmann
using ..CompositeGrids
using ..BrillouinZoneMeshes

locate(m::AbstractGrid, pos) = CompositeGrids.Interp.locate(m, pos)
volume(m::AbstractGrid, index) = CompositeGrids.Interp.volume(m, index)
volume(m::AbstractGrid) = CompositeGrids.Interp.volume(m)

locate(m::BrillouinZoneMeshes.AbstractMesh, pos) = BrillouinZoneMeshes.AbstractMeshes.locate(m, pos)
volume(m::BrillouinZoneMeshes.AbstractMesh, index) = BrillouinZoneMeshes.AbstractMeshes.volume(m, index)
volume(m::BrillouinZoneMeshes.AbstractMesh) = BrillouinZoneMeshes.AbstractMeshes.volume(m)

export locate, volume

abstract type TemporalGrid{T,REV} <: AbstractGrid{T} end

const FERMION = true
const BOSON = false
export FERMION, BOSON

export TemporalGrid
include("common.jl")

include("dlrfreq.jl")
export DLRFreq

include("imtime.jl")
export ImTime

include("imfreq.jl")
export ImFreq
export int_to_matfreq, matfreq_to_int
export matfreq

include("MeshProduct.jl")
export MeshProduct

#TODO: more functions from CompositeGrids


end
