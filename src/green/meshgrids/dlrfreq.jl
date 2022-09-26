"""
    struct DLRFreq{T, Grid} <: TemporalGrid{Int}

Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# Parameters
- `Grid`: type of 1D grid with T as the grid point type

# Members
- `grid`: 1D grid of time axis, with locate, volume, 
and AbstractArray interface implemented.
grid should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
"""
struct DLRFreq{T<:Real} <: TemporalGrid{T}
    dlr::DLRGrid
    grid::SimpleG.Arbitrary{T}
    beta::T
    Euv::T
    rtol::Float64
    sym::Symbol
    statistics::Statistics
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
function DLRFreq(beta, statistics::Statistics=UNKNOWN;
    dtype=Float64,
    rtol=1e-12,
    Euv=1000 / beta,
    sym=:none,
    dlr::Union{DLRGrid,Nothing}=nothing
)
    if isnothing(dlr)
        dlr = DLRGrid(Euv, beta, rtol, statistics isa Fermi, sym)
    end
    grid = SimpleG.Arbitrary{dtype}(dlr.ω)
    return DLRFreq{dtype}(dlr, grid, beta, Euv, rtol, sym, statistics)
end

Base.show(io::IO, tg::DLRFreq) = print(io, "DLR frequency grid with $(length(tg)) points, inverse temperature = $(tg.beta), UV Energy scale = $(tg.Euv), rtol = $(tg.rtol), sym = $(tg.sym), statistics = $(tg.statistics): $(_grid(tg.grid))")
