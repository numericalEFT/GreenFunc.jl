module MeshGrids

using ..GreenFunc
using ..GreenFunc: locate, volume
using ..Lehmann
using ..CompositeGrids

abstract type TemporalGrid{T} <: AbstractGrid{T} end

export TimeGrid

abstract type Statistics end
struct UnKnown <: Statistics end
struct Fermi <: Statistics end
struct Bose <: Statistics end
const UNKNOWN = UnKnown()
const FERMI = Fermi()
const BOSE = Bose()

Base.print(io::IO, s::UnKnown) = print(io, "UnKnown")
Base.print(io::IO, s::Fermi) = print(io, "Fermi")
Base.print(io::IO, s::Bose) = print(io, "Bose")
# Base.display(io::IO, s::Statistics) = Base.show(io, s)

# basic AbstractArray implement
Base.length(tg::TemporalGrid) = length(tg.grid)
Base.size(tg::TemporalGrid) = size(tg.grid)
Base.size(tg::TemporalGrid, I::Int) = size(tg.grid, I)
Base.getindex(tg::TemporalGrid, I::Int) = tg.grid[I]
Base.firstindex(tg::TemporalGrid) = 1
Base.lastindex(tg::TemporalGrid) = length(tg)

# iterator
Base.iterate(tg::TemporalGrid) = (tg[1], 1)
Base.iterate(tg::TemporalGrid, state) = (state >= length(tg)) ? nothing : (tg[state+1], state + 1)
# Base.IteratorSize(tg)
Base.IteratorSize(::Type{TemporalGrid{GT}}) where {GT} = Base.HasLength()
Base.IteratorEltype(::Type{TemporalGrid{GT}}) where {GT} = Base.HasEltype()
Base.eltype(::Type{TemporalGrid{GT}}) where {GT} = eltype(GT)

# locate and volume could fail if tg.grid has no implementation
GreenFunc.volume(tg::TemporalGrid, I::Int) = volume(tg.grid, I)
GreenFunc.volume(tg::TemporalGrid) = volume(tg.grid)
GreenFunc.locate(tg::TemporalGrid, pos) = locate(tg.grid, pos)

Base.floor(tg::TemporalGrid, pos) = floor(tg.grid, pos)

function _round(grid, sigdigits)
    if sigdigits <= 0
        return grid
    else
        return [round(x, sigdigits=sigdigits) for x in grid]
    end
end

# return a pretty print for grid. 
# If grid is very long, return [grid[1], grid[2], grid[3], ..., grid[end-2], grid[end-1], grid[end]]
# If grid is short, return the entire grid
function _grid(grid, n=3)
    if eltype(grid) <: Int
        # return join(io, ["[", join([grid[1:n]..., "...", grid[end-n:end]...], ", "), "]"])
        digits = 0
    else
        resolution = grid[2] - grid[1]
        digits = Int(round(log(grid[end] / resolution) / log(10))) + 3
        digits = digits < 5 ? 5 : digits
    end
    if length(grid) <= 2n + 3
        return join(["[", join(_round(grid, digits), ", "), "]"])
    else
        return join(["[", join([_round(grid[1:n], digits)..., "...", _round(grid[end-(n-1):end], digits)...], ", "), "]"])
    end
end

Base.show(io::IO, ::MIME"text/plain", tg::TemporalGrid) = Base.show(io, tg)
Base.show(io::IO, ::MIME"text/html", tg::TemporalGrid) = Base.show(io, tg)

#TODO: more functions from CompositeGrids

include("imtime.jl")
export ImTime


# meshgrids for green functions
# implement timegrids for time and frequency
# TODO: move MeshProduct into this module

"""
    struct TimeGrid{Domain<:TimeDomain,GridType}

Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# Parameters
- `Domain`: type of time domain, `Domain`<:`TimeDomain`.
- `GridType`: type of 1D grid

# Members
- `grid`: 1D grid of time axis, with locate, volume, 
and AbstractArray interface implemented.
grid should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
"""
struct TimeGrid{Domain<:TimeDomain,GridType}
    grid::GridType
    β::Float64
    isFermi::Bool
end

"""
    function TimeGrid(domain::Type{D<:TimeDomain};
         grid::GT, β, isFermi) where {D,GT}

Constructor of TimeGrid.

# Arguments
- `domain`: domain of time grid, Domain<:TimeDomain. By default, `domain` = IMFREQ.
- `grid`: time grid. All domain available for DLRGrid.
- `β`: Inverse temperature.
- `isFermi`: true if fermion, false if boson.
"""
function TimeGrid(::Type{D};
    grid::GT, β, isFermi) where {D,GT}
    return TimeGrid{D,GT}(grid, β, isFermi)
end

# for DLRGrid, default is DLRFreq with full DLRGrid stored
TimeGrid(grid::DLRGrid) = TimeGrid(
    DLRFreq, grid=grid, β=grid.β, isFermi=grid.isFermi)
# for ImTime and ImFreq, store only the corresponding grid
TimeGrid(::Type{<:ImTime}, grid::DLRGrid) = TimeGrid(
    ImTime, grid=CompositeGrids.SimpleG.Arbitrary{Float64}(grid.τ),
    β=grid.β, isFermi=grid.isFermi)
TimeGrid(::Type{<:ImFreq}, grid::DLRGrid) = TimeGrid(
    ImFreq, grid=CompositeGrids.SimpleG.Arbitrary{Float64}(grid.n),
    β=grid.β, isFermi=grid.isFermi)

# basic AbstractArray implement
Base.length(tg::TimeGrid) = length(tg.grid)
Base.size(tg::TimeGrid) = size(tg.grid)
Base.size(tg::TimeGrid, I::Int) = size(tg.grid, I)
Base.getindex(tg::TimeGrid, I::Int) = tg.grid[I]
Base.getindex(tg::TimeGrid{DLRFreq,DLRGrid}, I::Int) = tg.grid.ω[I]
Base.firstindex(tg::TimeGrid) = 1
Base.lastindex(tg::TimeGrid) = length(tg)

# iterator
Base.iterate(tg::TimeGrid) = (tg[1], 1)
Base.iterate(tg::TimeGrid, state) = (state >= length(tg)) ? nothing : (tg[state+1], state + 1)
# Base.IteratorSize(tg)
Base.IteratorSize(::Type{TimeGrid{D,GT}}) where {D,GT} = Base.HasLength()
Base.IteratorEltype(::Type{TimeGrid{D,GT}}) where {D,GT} = Base.HasEltype()
Base.eltype(::Type{TimeGrid{D,GT}}) where {D,GT} = eltype(GT)
Base.eltype(::Type{TimeGrid{D,DLRGrid}}) where {D} = eltype(fieldtypes(DLRGrid)[8])

ωn(tg::TimeGrid{ImFreq,GT}, I::Int) where {GT} = (tg.isFermi) ? π / tg.β * (2 * tg[I] + 1) : 2 * π / tg.β * tg[I]

# locate and volume could fail if tg.grid has no implementation
GreenFunc.volume(tg::TimeGrid, I::Int) = volume(tg.grid, I)
GreenFunc.volume(tg::TimeGrid) = volume(tg.grid)
GreenFunc.locate(tg::TimeGrid, pos) = locate(tg.grid, pos)

end