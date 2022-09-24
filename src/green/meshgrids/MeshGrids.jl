module MeshGrids

using ..GreenFunc
using ..Lehmann
using ..CompositeGrids

locate(m::AbstractGrid, pos) = CompositeGrids.Interp.locate(m, pos)
volume(m::AbstractGrid, index) = CompositeGrids.Interp.volume(m, index)
volume(m::AbstractGrid) = CompositeGrids.Interp.volume(m)

export locate, volume

abstract type TemporalGrid{T} <: AbstractGrid{T} end

abstract type Statistics end
struct UnKnown <: Statistics end
struct Fermi <: Statistics end
struct Bose <: Statistics end
const UNKNOWN = UnKnown()
const FERMI = Fermi()
const BOSE = Bose()

export UNKNOWN, FERMI, BOSE

include("common.jl")

include("imtime.jl")
# export ImTime

include("imfreq.jl")
# export ImFreq

include("dlrfreq.jl")
export DLRFreq

include("MeshProduct.jl")
export MeshProduct

#TODO: more functions from CompositeGrids

# export TimeGrid

# meshgrids for green functions
# implement timegrids for time and frequency
# # TODO: move MeshProduct into this module

# """
#     struct TimeGrid{Domain<:TimeDomain,GridType}

# Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# # Parameters
# - `Domain`: type of time domain, `Domain`<:`TimeDomain`.
# - `GridType`: type of 1D grid

# # Members
# - `grid`: 1D grid of time axis, with locate, volume, 
# and AbstractArray interface implemented.
# grid should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
# """
# struct TimeGrid{Domain<:TimeDomain,GridType}
#     grid::GridType
#     β::Float64
#     isFermi::Bool
# end

# """
#     function TimeGrid(domain::Type{D<:TimeDomain};
#          grid::GT, β, isFermi) where {D,GT}

# Constructor of TimeGrid.

# # Arguments
# - `domain`: domain of time grid, Domain<:TimeDomain. By default, `domain` = IMFREQ.
# - `grid`: time grid. All domain available for DLRGrid.
# - `β`: Inverse temperature.
# - `isFermi`: true if fermion, false if boson.
# """
# function TimeGrid(::Type{D};
#     grid::GT, β, isFermi) where {D,GT}
#     return TimeGrid{D,GT}(grid, β, isFermi)
# end

# # for DLRGrid, default is DLRFreq with full DLRGrid stored
# TimeGrid(grid::DLRGrid) = TimeGrid(
#     DLRFreq, grid=grid, β=grid.β, isFermi=grid.isFermi)
# # for ImTime and ImFreq, store only the corresponding grid
# TimeGrid(::Type{<:ImTime}, grid::DLRGrid) = TimeGrid(
#     ImTime, grid=CompositeGrids.SimpleG.Arbitrary{Float64}(grid.τ),
#     β=grid.β, isFermi=grid.isFermi)
# TimeGrid(::Type{<:ImFreq}, grid::DLRGrid) = TimeGrid(
#     ImFreq, grid=CompositeGrids.SimpleG.Arbitrary{Float64}(grid.n),
#     β=grid.β, isFermi=grid.isFermi)

# # basic AbstractArray implement
# Base.length(tg::TimeGrid) = length(tg.grid)
# Base.size(tg::TimeGrid) = size(tg.grid)
# Base.size(tg::TimeGrid, I::Int) = size(tg.grid, I)
# Base.getindex(tg::TimeGrid, I::Int) = tg.grid[I]
# Base.getindex(tg::TimeGrid{DLRFreq,DLRGrid}, I::Int) = tg.grid.ω[I]
# Base.firstindex(tg::TimeGrid) = 1
# Base.lastindex(tg::TimeGrid) = length(tg)

# # iterator
# Base.iterate(tg::TimeGrid) = (tg[1], 1)
# Base.iterate(tg::TimeGrid, state) = (state >= length(tg)) ? nothing : (tg[state+1], state + 1)
# # Base.IteratorSize(tg)
# Base.IteratorSize(::Type{TimeGrid{D,GT}}) where {D,GT} = Base.HasLength()
# Base.IteratorEltype(::Type{TimeGrid{D,GT}}) where {D,GT} = Base.HasEltype()
# Base.eltype(::Type{TimeGrid{D,GT}}) where {D,GT} = eltype(GT)
# Base.eltype(::Type{TimeGrid{D,DLRGrid}}) where {D} = eltype(fieldtypes(DLRGrid)[8])

# ωn(tg::TimeGrid{ImFreq,GT}, I::Int) where {GT} = (tg.isFermi) ? π / tg.β * (2 * tg[I] + 1) : 2 * π / tg.β * tg[I]

# locate and volume could fail if tg.grid has no implementation
# GreenFunc.volume(tg::TimeGrid, I::Int) = volume(tg.grid, I)
# GreenFunc.volume(tg::TimeGrid) = volume(tg.grid)
# GreenFunc.locate(tg::TimeGrid, pos) = locate(tg.grid, pos)

end