"""
    struct ImFreq{T, Grid} <: TemporalGrid{Int}

Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# Parameters
- `Grid`: type of 1D grid with T as the grid point type

# Members
- `grid`: 1D grid of time axis, with locate, volume, 
and AbstractArray interface implemented.
grid should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
"""
struct ImFreq{T<:Real,Grid} <: TemporalGrid{Int}
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
function ImFreq(beta, statistics::Statistics=UNKNOWN;
             dtype=Int64,
    rtol = 1e-12,
    Euv=1000 / beta,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    if isnothing(grid)
        dlr = DLRGrid(Euv, beta, rtol, statistics isa Fermi, :none)
        grid = SimpleG.Arbitrary{Int}(dlr.n)
    elseif (grid isa AbstractVector)
        grid = SimpleG.Arbitrary{Int}(Int.(grid))
    end
    return ImFreq{dtype,typeof(grid)}(grid, beta, Euv, statistics)
end

matfreq_to_int(tg::ImFreq, ωn) = (tg.statistics isa Fermi) ? Int(round((ωn * tg.beta / π - 1) / 2)) : Int(round((ωn * tg.beta / π) / 2))
int_to_matfreq(tg::ImFreq, n::Int) = (tg.statistics isa Fermi) ? (2n + 1) * π / tg.beta : 2n * π / tg.beta

"""
    Base.getindex(g::ImFreq, I::Int)

Equivalent to `g[I]`, get the __real-valued__ Matsubara frequency of the Ith point in the grid. 
For fermion, return (2g[I]+1)π/β, for boson, return 2g[I]*π/β.

If you need the __integer-valued__ frequency, use `g.grid[I]` instead.
"""
Base.getindex(tg::ImFreq, I::Int) = int_to_matfreq(tg, tg.grid[I])
Base.show(io::IO, tg::ImFreq) = print(io, "Matsubara frequency grid with $(length(tg)) points, inverse temperature = $(tg.beta), UV Energy scale = $(tg.Euv), statistics = $(tg.statistics): $(_grid(tg.grid))")


locate(tg::ImFreq, n::Int) = locate(tg.grid, n)

function locate(tg::ImFreq, ωn)
    n = matfreq_to_int(tg, ωn)
    return locate(tg.grid, n)
end
