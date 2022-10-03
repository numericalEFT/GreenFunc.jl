module MeshGrids

using ..GreenFunc
using ..Lehmann
using ..CompositeGrids
using ..BrillouinZoneMeshes

locate(m::AbstractGrid, pos) = CompositeGrids.Interp.locate(m, pos)
volume(m::AbstractGrid, index) = CompositeGrids.Interp.volume(m, index)
volume(m::AbstractGrid) = CompositeGrids.Interp.volume(m)

locate(m::BrillouinZoneMeshes.BaseMesh.AbstractMesh, pos) = BrillouinZoneMeshes.BaseMesh.locate(m, pos)
volume(m::BrillouinZoneMeshes.BaseMesh.AbstractMesh, index) = BrillouinZoneMeshes.BaseMesh.volume(m, index)
volume(m::BrillouinZoneMeshes.BaseMesh.AbstractMesh) = BrillouinZoneMeshes.BaseMesh.volume(m, pos)



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
