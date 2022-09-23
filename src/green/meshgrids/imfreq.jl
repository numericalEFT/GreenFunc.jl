"""
    struct ImFreq{Grid} <: TemporalGrid{Int}

Time grid for Green's functions. An 1D grid wrapped with TimeDomain tpye.

# Parameters
- `Grid{T}`: type of 1D grid with T as the grid point type

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
    dtype=Float64,
    Euv=128 / beta,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    if isnothing(grid)
        dlr = DLRGrid(Euv, beta, rtol, statistics isa Fermi, :none)
        println(dlr)
        grid = SimpleG.Arbitrary{Int}(dlr.n)
    elseif (grid isa AbstractVector)
        grid = SimpleG.Arbitrary{Int}(Int.(grid))
    end
    return ImFreq{dtype,typeof(grid)}(grid, beta, Euv, statistics)
end

Base.getindex(tg::ImFreq, I::Int) = (tg.statistics isa Fermi) ? π / tg.beta * (2 * tg.grid[I] + 1) : 2 * π / tg.beta * tg.grid[I]

Base.show(io::IO, tg::ImFreq) = print(io, "Matsubara frequency grid with $(length(tg)) points, inverse temperature = $(tg.beta), UV Energy scale = $(tg.Euv), statistics = $(tg.statistics): $(_grid(tg.grid))")


# ωn(tg::ImFreq{GT}, I::Int) where {GT} = (tg.statistics isa Fermi) ? π / tg.beta * (2 * tg[I] + 1) : 2 * π / tg.beta * tg[I]
# matfreq(tg::ImFreq, I::Int) = ωn(tg, I)