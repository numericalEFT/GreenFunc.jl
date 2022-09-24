"""
    struct ImTime{T, Grid} <: TemporalGrid{T}

Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# Parameters
- `T`: type of the grid point
- `Grid`: type of 1D grid with T as the grid point type

# Members
- `grid`: 1D grid of time axis, with locate, volume, 
and AbstractArray interface implemented.
grid should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
"""
struct ImTime{T<:Real,Grid} <: TemporalGrid{T}
    grid::Grid
    beta::T
    Euv::T
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
function ImTime(beta, statistics::Statistics=UNKNOWN;
    dtype=Float64,
    rtol=1e-12,
    Euv=1000 / beta,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    if isnothing(grid)
        dlr = DLRGrid(Euv, beta, rtol, statistics isa Fermi, :none)
        grid = SimpleG.Arbitrary{dtype}(dlr.τ)
        # grid = SimpleG.Uniform{dtype}([0, beta], Int(round(beta / resolution)))
        # grid = CompositeGrid.LogDensedGrid(:uniform, [0.0, beta], [0.0, beta], 8, 1 / Euv, 8) #roughly ~100 points if resolution = beta/128
    elseif (grid isa AbstractVector)
        grid = SimpleG.Arbitrary{dtype}(grid)
    end
    @assert eltype(grid) == dtype "The type of grid should be the same as dtype = $dtype"
    return ImTime{dtype,typeof(grid)}(grid, beta, Euv, statistics)
end

Base.show(io::IO, tg::ImTime) = print(io, "Imaginary Time grid with $(length(tg)) points, inverse temperature = $(tg.beta), UV Energy scale = $(tg.Euv), statistics = $(tg.statistics): $(_grid(tg.grid))")
