"""
    struct DLRFreq{T<:Real} <: TemporalGrid{Int}

Discrete-Lehmann-representation grid for Green's functions. 

# Parameters
- `T`: type of the `grid` point, `β` and `Euv`.

# Members
- `dlr`: built-in DLR grid.
- `grid`: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented.
  It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
- `β`: inverse temperature.
- `Euv`:  the UV energy scale of the spectral density.
- `rtol`: tolerance absolute error.
- `sym`: the symmetry of `dlr`.
- `statistics`: type of statistics for particles. It can be `FERMI`, `BOSE`, and `UNKNOWN`.
"""
struct DLRFreq{T<:Real} <: TemporalGrid{T}
    dlr::DLRGrid
    grid::SimpleG.Arbitrary{T}
    β::T
    Euv::T
    rtol::Float64
    sym::Symbol
    statistics::Statistics
end

"""
    function DLRFreq(β, statistics::Statistics=UNKNOWN;
        dtype=Float64,
        rtol=1e-12,
        Euv=1000 / β,
        sym=:none,
        dlr::Union{DLRGrid,Nothing}=nothing
    )

Create a `DLRFreq` struct.

# Arguments
- `β`: inverse temperature.
- `statistics`: type of statistics for particles, including `FERMI`, `BOSE`, and `UNKNOWN`. By default, `statistics = UNKNOWN`.
- `dtype`: type of `β` and `Euv`.
- `rtol`: tolerance absolute error. By default, `rtol` = 1e-12.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv = 1000 / β`.
- `sym`: the symmetry of `dlr`. By default, `sym = :none`.
- `dlr`: 1D DLR grid. By default, a DLR grid with input arguments is used.
"""
function DLRFreq(β, statistics::Statistics=UNKNOWN;
    dtype=Float64,
    rtol=1e-12,
    Euv=1000 / β,
    sym=:none,
    dlr::Union{DLRGrid,Nothing}=nothing
)
    if isnothing(dlr)
        dlr = DLRGrid(Euv, β, rtol, statistics isa Fermi, sym)
    end
    grid = SimpleG.Arbitrary{dtype}(dlr.ω)
    return DLRFreq{dtype}(dlr, grid, β, Euv, rtol, sym, statistics)
end

"""
    show(io::IO, tg::DLRFreq)

Write a text representation of the DLR grid `tg` to the output stream `io`.
"""
Base.show(io::IO, tg::DLRFreq) = print(io, "DLR frequency grid with $(length(tg)) points, inverse temperature = $(tg.β), UV Energy scale = $(tg.Euv), rtol = $(tg.rtol), sym = $(tg.sym), statistics = $(tg.statistics): $(_grid(tg.grid))")